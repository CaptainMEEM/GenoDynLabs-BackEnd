"""
The unit of work, isolated from any web framework or queue. Call it from an
RQ/Celery worker, a thread, or directly. Heavy reference data is loaded once
per process (inside annotate) and reused across every job.
"""
import os

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

BASE_DIR  = os.path.dirname(os.path.abspath(__file__))
DATA_DIR  = os.path.join(BASE_DIR, "data")
REF_PKL   = os.path.join(DATA_DIR, "reference.pkl")
TRAIT_CSV = os.path.join(DATA_DIR, "trait_df.csv")
EQ_CSV    = os.path.join(DATA_DIR, "equilibrium_df.csv")


def run_report_job(job_dir, clean_path, user_email, display_name=""):
    """Annotate -> build PDF -> email. Returns variant count. Safe to retry."""
    snpedia_csv = os.path.join(job_dir, "snpedia_data.csv")
    n = annotate_genome(clean_path, snpedia_csv, ref_path=REF_PKL)   # ~0.3s

    pdf = build_pdf(snpedia_csv, trait_csv=TRAIT_CSV, eq_csv=EQ_CSV,
                    user_display_name=display_name, require_note=True)
    with open(os.path.join(job_dir, "genodynlabs_report.pdf"), "wb") as f:
        f.write(pdf)

    if send_report_email and user_email:
        send_report_email(to_address=user_email, user_display_name=display_name,
                          pdf_bytes=pdf, pdf_filename="genodynlabs_report.pdf")
    return n
