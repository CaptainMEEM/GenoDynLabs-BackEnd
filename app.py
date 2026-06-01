import os
import sys
import subprocess
import uuid
import logging
from flask import Flask, request, jsonify, send_file
from flask_cors import CORS

from services.firebase_auth import init_firebase, require_auth
from services.genodyn_report import build_pdf          # condensed, topic-organized report
from services.email_sender  import send_report_email

# ── Setup ──────────────────────────────────────────────────────────────
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("app")

app = Flask(__name__)

# Restrict CORS to your Vercel frontend in production.
# Set ALLOWED_ORIGINS in Railway, e.g. "https://genodynlabs.vercel.app,https://genodynlabs.com"
ALLOWED_ORIGINS = [o.strip() for o in os.environ.get("ALLOWED_ORIGINS", "*").split(",")]
CORS(app, resources={r"/*": {"origins": ALLOWED_ORIGINS}}, supports_credentials=False)

BASE_DIR    = os.path.dirname(os.path.abspath(__file__))
UPLOAD_DIR  = os.path.join(BASE_DIR, "uploads")
OUTPUT_DIR  = os.path.join(BASE_DIR, "outputs")
SCRIPT      = os.path.join(BASE_DIR, "generate_variant_report.py")
TRAIT_CSV   = os.path.join(BASE_DIR, "data", "trait_df.csv")        # GWAS effect alleles + odds ratios
EQ_CSV      = os.path.join(BASE_DIR, "data", "equilibrium_df.csv")  # LD r^2 for dedup
ALLOWED_EXT = {".txt"}
MAX_MB      = 50

app.config["MAX_CONTENT_LENGTH"] = MAX_MB * 1024 * 1024
os.makedirs(UPLOAD_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Initialize Firebase on startup so misconfiguration fails fast.
try:
    init_firebase()
    log.info("Firebase Admin SDK initialized")
except Exception as e:
    log.warning(f"Firebase init deferred (will retry on first request): {e}")


def allowed(filename):
    return os.path.splitext(filename)[1].lower() in ALLOWED_EXT


def preprocess_23andme(src_path, dst_path):
    """Strip 23andMe comment lines, normalize header."""
    with open(src_path, "r", errors="replace") as f:
        lines = f.readlines()

    out = []
    header_written = False
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            if "rsid" in stripped and not header_written:
                out.append("rsid\tchromosome\tposition\tgenotype\n")
                header_written = True
        else:
            if not header_written:
                out.append(line)
                header_written = True
            else:
                out.append(line)

    with open(dst_path, "w") as f:
        f.writelines(out)


# ── Routes ─────────────────────────────────────────────────────────────

@app.route("/health")
def health():
    return jsonify({"status": "ok"})


@app.route("/upload", methods=["POST"])
@require_auth
def upload(user):
    """
    Verified Firebase user uploads their raw 23andMe .txt file.
    We run Modern Promethease, build a PDF, and email it to their
    verified Firebase email.
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

    job_id  = str(uuid.uuid4())
    job_dir = os.path.join(OUTPUT_DIR, job_id)
    os.makedirs(job_dir, exist_ok=True)
    log.info(f"[{job_id}] upload from uid={uid} email={user_email}")

    raw_path   = os.path.join(UPLOAD_DIR, f"{job_id}_raw.txt")
    clean_path = os.path.join(job_dir, "genome_clean.txt")
    file.save(raw_path)

    try:
        preprocess_23andme(raw_path, clean_path)

        with open(clean_path) as f:
            header = f.readline().strip().lower()
        if "rsid" not in header:
            return jsonify({"error": "File does not look like a valid 23andMe file"}), 400

        # Run Modern Promethease
        log.info(f"[{job_id}] running Modern Promethease")
        result = subprocess.run(
            [sys.executable, SCRIPT, clean_path],
            capture_output=True, text=True, cwd=job_dir, timeout=300,
        )
        if result.returncode != 0:
            log.error(f"[{job_id}] promethease failed: {result.stderr[-500:]}")
            return jsonify({
                "error": "Promethease processing failed",
                "details": result.stderr[-2000:],
            }), 500

        snpedia_csv = os.path.join(job_dir, "snpedia_data.csv")
        if not os.path.exists(snpedia_csv):
            return jsonify({"error": "Annotated CSV was not generated"}), 500

        # Build PDF (condensed, topic-organized; enriched with GWAS traits + LD dedup)
        log.info(f"[{job_id}] building PDF")
        pdf_bytes = build_pdf(
            snpedia_csv,
            trait_csv=TRAIT_CSV,
            eq_csv=EQ_CSV,
            user_display_name=display_name,
            drop_neutral_zero=True,
        )
        pdf_path = os.path.join(job_dir, "genodynlabs_report.pdf")
        with open(pdf_path, "wb") as f:
            f.write(pdf_bytes)

        # Email it
        log.info(f"[{job_id}] emailing PDF to {user_email}")
        send_report_email(
            to_address=user_email,
            user_display_name=display_name,
            pdf_bytes=pdf_bytes,
            pdf_filename="genodynlabs_report.pdf",
        )

        # Count variants for the response
        with open(snpedia_csv) as f:
            total = sum(1 for _ in f) - 1

        return jsonify({
            "job_id":         job_id,
            "total_variants": total,
            "emailed_to":     user_email,
            "message":        "Your report has been emailed to you.",
        })

    except subprocess.TimeoutExpired:
        log.error(f"[{job_id}] promethease timed out")
        return jsonify({"error": "Processing timed out"}), 500
    except Exception as e:
        log.exception(f"[{job_id}] failed")
        return jsonify({"error": str(e)}), 500
    finally:
        if os.path.exists(raw_path):
            os.remove(raw_path)


if __name__ == "__main__":
    app.run(debug=True, port=5000)
