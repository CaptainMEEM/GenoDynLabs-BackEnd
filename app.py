import os
import uuid
import logging
import threading
from concurrent.futures import ThreadPoolExecutor

from flask import Flask, request, jsonify
from flask_cors import CORS

from services.firebase_auth import init_firebase, require_auth
from services.jobs import run_report_job

# ── Setup ──────────────────────────────────────────────────────────────
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("app")

app = Flask(__name__)

# Restrict CORS to your Vercel frontend in production.
# Set ALLOWED_ORIGINS in Railway, e.g. "https://genodynlabs.vercel.app,https://genodynlabs.com"
ALLOWED_ORIGINS = [o.strip() for o in os.environ.get("ALLOWED_ORIGINS", "*").split(",")]
CORS(app, resources={r"/*": {"origins": ALLOWED_ORIGINS}}, supports_credentials=False)

ALLOWED_EXT = {".txt"}
MAX_MB      = 50
app.config["MAX_CONTENT_LENGTH"] = MAX_MB * 1024 * 1024

# Initialize Firebase on startup so misconfiguration fails fast.
try:
    init_firebase()
    log.info("Firebase Admin SDK initialized")
except Exception as e:
    log.warning(f"Firebase init deferred (will retry on first request): {e}")

# ── Job dispatch: RQ if Redis is configured, else a local thread pool ────
# The caller (/upload) uses the same dispatch() either way; only the executor
# differs. Add the Redis plugin + set REDIS_URL and this upgrades to the real
# worker (worker.py) with no further code change. Until then jobs run in a
# background thread so /upload still returns immediately.
REDIS_URL   = os.environ.get("REDIS_URL")
QUEUE_NAME  = "reports"     # must match worker.py
JOB_TIMEOUT = 600           # seconds; annotate is fast, PDF + email dominate
RESULT_TTL  = 3600          # keep finished/failed records this long for polling

if REDIS_URL:
    import redis
    from rq import Queue
    from rq.job import Job
    from rq.exceptions import NoSuchJobError

    _redis = redis.from_url(REDIS_URL)
    _queue = Queue(QUEUE_NAME, connection=_redis)
    log.info("dispatch: RQ via REDIS_URL")
else:
    # Fallback for local dev / no-Redis deploys. NOTE: status lives in this
    # process's memory, so it is reliable only with a single web process
    # (the dev server, or gunicorn --workers 1). The email is still sent from
    # whichever process ran the job, so delivery is unaffected by worker count.
    _executor      = ThreadPoolExecutor(max_workers=2)
    _jobs_lock     = threading.Lock()
    _jobs          = {}        # job_id -> {"status","uid","variants","error"}
    MAX_LOCAL_JOBS = 2000      # bound the in-memory store
    log.info("dispatch: in-process ThreadPoolExecutor (no REDIS_URL)")


def _local_runner(job_id, uid, genome_text, user_email, display_name):
    with _jobs_lock:
        _jobs[job_id] = {"status": "processing", "uid": uid}
    try:
        n = run_report_job(genome_text, user_email, display_name)
        with _jobs_lock:
            _jobs[job_id] = {"status": "done", "uid": uid, "variants": n}
        log.info(f"[{job_id}] done ({n} variants)")
    except Exception as e:
        log.exception(f"[{job_id}] failed")
        with _jobs_lock:
            _jobs[job_id] = {"status": "failed", "uid": uid, "error": str(e)}


def dispatch(job_id, uid, genome_text, user_email, display_name):
    if REDIS_URL:
        _queue.enqueue(
            run_report_job, genome_text, user_email, display_name,
            job_id=job_id, job_timeout=JOB_TIMEOUT,
            result_ttl=RESULT_TTL, failure_ttl=RESULT_TTL,
            meta={"uid": uid},
        )
    else:
        with _jobs_lock:
            if len(_jobs) >= MAX_LOCAL_JOBS:                 # prune oldest half
                for k in list(_jobs)[: len(_jobs) // 2]:
                    _jobs.pop(k, None)
            _jobs[job_id] = {"status": "queued", "uid": uid}
        _executor.submit(_local_runner, job_id, uid, genome_text,
                         user_email, display_name)


# ── 23andMe preprocessing (in memory, no temp files) ─────────────────────
def preprocess_23andme_text(raw_text):
    """Strip 23andMe comment lines and normalize the header. Pure string work."""
    out = []
    header_written = False
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


def allowed(filename):
    return os.path.splitext(filename)[1].lower() in ALLOWED_EXT


# ── Status vocabulary: queued | processing | done | failed | unknown ─────
# Same words for both backends so the frontend handles one set.
_RQ_STATUS = {
    "queued": "queued", "deferred": "queued", "scheduled": "queued",
    "started": "processing",
    "finished": "done",
    "failed": "failed", "stopped": "failed", "canceled": "failed",
}


# ── Routes ───────────────────────────────────────────────────────────────
@app.route("/health")
def health():
    return jsonify({"status": "ok", "dispatch": "rq" if REDIS_URL else "thread"})


@app.route("/upload", methods=["POST"])
@require_auth
def upload(user):
    """
    Cheap path only: auth, validate the file, strip 23andMe comments in memory,
    hand the job off, and return 202 + job_id immediately. Annotation, PDF
    rendering and email all happen in the worker. The email is the primary
    delivery; /status is for progress UI.
    """
    uid          = user.get("uid")
    user_email   = user.get("email")
    display_name = user.get("name") or ""

    if not user_email:
        return jsonify({"error": "Your account has no email on file"}), 400
    if "file" not in request.files:
        return jsonify({"error": "No file part in request"}), 400

    file = request.files["file"]
    if file.filename == "":
        return jsonify({"error": "No file selected"}), 400
    if not allowed(file.filename):
        return jsonify({"error": "Only .txt files are accepted"}), 400

    raw_text   = file.read().decode("utf-8", errors="replace")
    clean_text = preprocess_23andme_text(raw_text)

    first_line = clean_text.split("\n", 1)[0].lower()
    if "rsid" not in first_line:
        return jsonify({"error": "File does not look like a valid 23andMe file"}), 400
    if clean_text.count("\n") < 2:                     # header only, no data rows
        return jsonify({"error": "File contains no genotype rows"}), 400

    job_id = str(uuid.uuid4())
    log.info(f"[{job_id}] upload from uid={uid} email={user_email}; dispatching")
    dispatch(job_id, uid, clean_text, user_email, display_name)

    return jsonify({
        "job_id":     job_id,
        "status":     "queued",
        "emailed_to": user_email,
        "message":    "Your report is being generated and will be emailed to you.",
    }), 202


@app.route("/status/<job_id>")
@require_auth
def status(user, job_id):
    """
    Poll a job. Mirrors upload's auth wiring: `user` is injected by
    require_auth and `job_id` comes from the URL. Only the job's owner can
    read it (uid check) so status can't be probed across accounts.
    """
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
