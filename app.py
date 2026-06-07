"""
app.py  --  GenoDynLabs web service (cheap, non-blocking path ONLY).

The web process does the minimum and returns in milliseconds:
  1. authenticate (Firebase),
  2. read + validate the upload,
  3. (if sealed) open the client's sealed-box; else strip 23andMe comments,
  4. ENCRYPT the genome at rest,
  5. enqueue the job and return 202 + job_id.

All heavy work (annotate, render, email) happens in the RQ worker. With Redis
configured the web service holds no per-job state, so you can run multiple web
replicas behind Railway's load balancer and scale workers independently. A
thread-pool fallback exists only for single-process local dev.
"""
import os
import uuid
import logging

from flask import Flask, request, jsonify
from flask_cors import CORS

from services.firebase_auth import init_firebase, require_auth
from services.jobs import run_report_job
from services import crypto

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("app")

app = Flask(__name__)

ALLOWED_ORIGINS = [o.strip() for o in os.environ.get("ALLOWED_ORIGINS", "*").split(",")]
CORS(app, resources={r"/*": {"origins": ALLOWED_ORIGINS}}, supports_credentials=False)

ALLOWED_EXT = {".txt"}
MAX_MB = 60
app.config["MAX_CONTENT_LENGTH"] = MAX_MB * 1024 * 1024

try:
    init_firebase()
    log.info("Firebase Admin SDK initialized")
except Exception as e:
    log.warning(f"Firebase init deferred: {e}")

if not crypto.encryption_enabled():
    log.warning("GENODYN_DATA_KEYS not set — payloads will NOT be encrypted at "
                "rest. Set it in production.")

# ── Dispatch: RQ when REDIS_URL is set, else a local thread pool (dev only) ──
REDIS_URL = os.environ.get("REDIS_URL")
QUEUE_NAME = "reports"
JOB_TIMEOUT = 900
RESULT_TTL = 3600

if REDIS_URL:
    import redis
    from rq import Queue
    from rq.job import Job
    from rq.exceptions import NoSuchJobError
    _redis = redis.from_url(REDIS_URL)
    _queue = Queue(QUEUE_NAME, connection=_redis)
    log.info("dispatch: RQ via REDIS_URL")
else:
    import threading
    from concurrent.futures import ThreadPoolExecutor
    _executor = ThreadPoolExecutor(max_workers=2)
    _jobs_lock = threading.Lock()
    _jobs = {}
    MAX_LOCAL_JOBS = 2000
    log.info("dispatch: in-process ThreadPoolExecutor (no REDIS_URL)")

# Fail CLOSED on a real (queued) deployment: if the genome would be written to
# Redis, it MUST be encrypted at rest first. We refuse rather than silently
# persisting a plaintext genome. The local thread-pool dev path (no Redis, no
# persistence) is allowed to run unencrypted with the warning above.
REQUIRE_ENCRYPTION = bool(REDIS_URL)
if REQUIRE_ENCRYPTION and not crypto.encryption_enabled():
    log.error("REDIS_URL is set but GENODYN_DATA_KEYS is not — uploads will be "
              "refused until at-rest encryption is configured.")


def _local_runner(job_id, uid, enc_payload, user_email, display_name):
    with _jobs_lock:
        _jobs[job_id] = {"status": "processing", "uid": uid}
    try:
        n = run_report_job(enc_payload, user_email, display_name)
        with _jobs_lock:
            _jobs[job_id] = {"status": "done", "uid": uid, "variants": n}
    except Exception as e:
        log.exception(f"[{job_id}] failed")
        with _jobs_lock:
            _jobs[job_id] = {"status": "failed", "uid": uid, "error": str(e)}


def dispatch(job_id, uid, enc_payload, user_email, display_name):
    if REDIS_URL:
        _queue.enqueue(run_report_job, enc_payload, user_email, display_name,
                       job_id=job_id, job_timeout=JOB_TIMEOUT,
                       result_ttl=RESULT_TTL, failure_ttl=RESULT_TTL,
                       meta={"uid": uid})
    else:
        with _jobs_lock:
            if len(_jobs) >= MAX_LOCAL_JOBS:
                for k in list(_jobs)[: len(_jobs) // 2]:
                    _jobs.pop(k, None)
            _jobs[job_id] = {"status": "queued", "uid": uid}
        _executor.submit(_local_runner, job_id, uid, enc_payload,
                         user_email, display_name)


# ── 23andMe preprocessing (pure string work, in memory) ──────────────────
def preprocess_23andme_text(raw_text):
    out, header_written = [], False
    for line in raw_text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            if "rsid" in stripped and not header_written:
                out.append("rsid\tchromosome\tposition\tgenotype")
                header_written = True
        else:
            out.append(line)
            header_written = True
    return ("\n".join(out) + "\n") if out else ""


def _allowed(filename):
    return os.path.splitext(filename)[1].lower() in ALLOWED_EXT


_RQ_STATUS = {
    "queued": "queued", "deferred": "queued", "scheduled": "queued",
    "started": "processing", "finished": "done",
    "failed": "failed", "stopped": "failed", "canceled": "failed",
}


# ── Routes ───────────────────────────────────────────────────────────────
@app.route("/health")
def health():
    return jsonify({"status": "ok",
                    "dispatch": "rq" if REDIS_URL else "thread",
                    "encryption": crypto.encryption_enabled()})


@app.route("/pubkey")
def pubkey():
    """X25519 public key for optional client-side sealed-box uploads.
    Empty string means the sealed-box path isn't configured on this deploy."""
    return jsonify({"x25519_public_key": crypto.server_public_key_b64()})


@app.route("/upload", methods=["POST"])
@require_auth
def upload(user):
    uid = user.get("uid")
    user_email = user.get("email")
    display_name = user.get("name") or ""

    # Never persist a plaintext genome to the queue.
    if REQUIRE_ENCRYPTION and not crypto.encryption_enabled():
        log.error("refusing upload: at-rest encryption not configured")
        return jsonify({"error": "Service is temporarily unavailable "
                                 "(secure storage not configured)."}), 503

    if not user_email:
        return jsonify({"error": "Your account has no email on file"}), 400
    if "file" not in request.files:
        return jsonify({"error": "No file part in request"}), 400

    file = request.files["file"]
    if file.filename == "":
        return jsonify({"error": "No file selected"}), 400

    raw = file.read()

    # Path A: client sealed the genome to our public key (opaque upload).
    if crypto.is_sealed(raw):
        try:
            genome_bytes = crypto.open_sealed(raw)
        except Exception as e:
            log.warning(f"sealed-box open failed: {e}")
            return jsonify({"error": "Could not open the encrypted upload"}), 400
        raw_text = genome_bytes.decode("utf-8", errors="replace")
    else:
        # Path B: plaintext .txt over TLS.
        if not _allowed(file.filename):
            return jsonify({"error": "Only .txt files are accepted"}), 400
        raw_text = raw.decode("utf-8", errors="replace")

    clean_text = preprocess_23andme_text(raw_text)
    first_line = clean_text.split("\n", 1)[0].lower()
    if "rsid" not in first_line:
        return jsonify({"error": "File does not look like a valid 23andMe file"}), 400
    if clean_text.count("\n") < 2:
        return jsonify({"error": "File contains no genotype rows"}), 400

    # Encrypt at rest BEFORE the genome enters the queue.
    enc_payload = crypto.encrypt(clean_text.encode("utf-8"))

    job_id = str(uuid.uuid4())
    log.info(f"[{job_id}] upload uid={uid}; sealed={crypto.is_sealed(raw)}; dispatching")
    dispatch(job_id, uid, enc_payload, user_email, display_name)

    return jsonify({
        "job_id": job_id,
        "status": "queued",
        "emailed_to": user_email,
        "message": "Your report is being generated and will be emailed to you.",
    }), 202


@app.route("/status/<job_id>")
@require_auth
def status(user, job_id):
    uid = user.get("uid")
    if REDIS_URL:
        try:
            job = Job.fetch(job_id, connection=_redis)
        except NoSuchJobError:
            return jsonify({"job_id": job_id, "status": "unknown"}), 404
        if job.meta.get("uid") != uid:
            return jsonify({"job_id": job_id, "status": "unknown"}), 404
        state = _RQ_STATUS.get(job.get_status(refresh=True), "unknown")
        body = {"job_id": job_id, "status": state}
        if state == "done":
            body["variants"] = job.result
        elif state == "failed":
            body["error"] = "Report generation failed"
        return jsonify(body)

    with _jobs_lock:
        info = _jobs.get(job_id)
    if not info or info.get("uid") != uid:
        return jsonify({"job_id": job_id, "status": "unknown"}), 404
    body = {"job_id": job_id, "status": info["status"]}
    if info["status"] == "done":
        body["variants"] = info.get("variants")
    elif info["status"] == "failed":
        body["error"] = "Report generation failed"
    return jsonify(body)


if __name__ == "__main__":
    app.run(debug=True, port=5000)
