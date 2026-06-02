"""
The unit of work, isolated from any web framework or queue. Call it from an
RQ/Celery worker, a thread, or directly. Heavy reference data is loaded once
per process (inside annotate) and reused across every job.

The job receives the *cleaned genome text*, not a path. On Railway the web and
worker run as separate services with separate filesystems, so a path written by
the web dyno is invisible to the worker. The text travels in the job payload
(through Redis). Each job writes its own scratch dir and removes it when done,
so nothing accumulates on the worker and retries stay safe.
"""
import os
import shutil
import tempfile

try:
    from . import annotate, genodyn_report
    from .annotate import annotate_genome
    from .genodyn_report import build_pdf
    from .email_sender import send_report_email
except ImportError:                      # standalone / testing
    from annotate import annotate_genome
    from genodyn_report import build_pdf
    try:
        from email_sender import send_report_email
    except ImportError:
        send_report_email = None

BASE_DIR = os.path.dirname(os.path.abspath(__file__))    # .../services


def _find_data_dir():
    """Locate the folder holding reference.pkl. Tries root-level data/ (this
    repo's layout) first, then services/data/, so it works regardless of where
    the data lands. Falls back to root data/ if the pkl isn't found anywhere."""
    candidates = [
        os.path.join(os.path.dirname(BASE_DIR), "data"),  # <repo root>/data   (your layout)
        os.path.join(BASE_DIR, "data"),                   # services/data
    ]
    for d in candidates:
        if os.path.exists(os.path.join(d, "reference.pkl")):
            return d
    return candidates[0]


DATA_DIR  = _find_data_dir()
REF_PKL   = os.path.join(DATA_DIR, "reference.pkl")
TRAIT_CSV = os.path.join(DATA_DIR, "trait_df.csv")
EQ_CSV    = os.path.join(DATA_DIR, "equilibrium_df.csv")


def run_report_job(genome_text, user_email, display_name=""):
    """Annotate -> build PDF -> email. Returns variant count. Safe to retry.

    `genome_text` is the comment-stripped, header-normalized 23andMe file as a
    single string (the web service already did that cheap part). We materialize
    it in a private scratch dir because annotate_genome / build_pdf operate on
    file paths, then tear the dir down regardless of outcome.
    """
    scratch = tempfile.mkdtemp(prefix="genodyn_")
    try:
        genome_path = os.path.join(scratch, "genome_clean.txt")
        with open(genome_path, "w", newline="") as f:
            f.write(genome_text)

        snpedia_csv = os.path.join(scratch, "snpedia_data.csv")
        n = annotate_genome(genome_path, snpedia_csv, ref_path=REF_PKL)   # ~0.3s

        pdf = build_pdf(snpedia_csv, trait_csv=TRAIT_CSV, eq_csv=EQ_CSV,
                        user_display_name=display_name, require_note=True)

        if send_report_email and user_email:
            send_report_email(to_address=user_email, user_display_name=display_name,
                              pdf_bytes=pdf, pdf_filename="genodynlabs_report.pdf")
        return n
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
