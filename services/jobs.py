"""
jobs.py  --  The unit of work, isolated from any web framework or queue.

Design goals this file delivers on:
  * NON-BLOCKING / mass scale -- this function is what the RQ worker runs off
    the queue; the web process never does any of this heavy work. Many workers
    can run it concurrently (it shares only the read-only, process-cached
    reference bundle).
  * FULLY IN MEMORY -- the genome is decrypted into memory, annotated, rendered,
    emailed, and discarded. Nothing (no genome, no PDF) is ever written to disk,
    so there is no plaintext at rest and nothing to clean up between jobs.
  * ENCRYPTED PAYLOAD -- the argument is the AES-256-GCM-wrapped genome produced
    by app.py, so the genome sits in Redis only as ciphertext.
  * LOW RAILWAY COST -- the reference bundle loads once per worker process and is
    reused for every job; per-job allocation is just the user's records.
"""
import os
import logging

try:
    from . import crypto
    from .annotate import annotate_genome, load_reference
    from .report import build_pdf
    from .email_sender import send_report_email
except ImportError:                      # standalone / testing
    import crypto
    from annotate import annotate_genome, load_reference
    from report import build_pdf
    try:
        from email_sender import send_report_email
    except ImportError:
        send_report_email = None

log = logging.getLogger("jobs")

BASE_DIR = os.path.dirname(os.path.abspath(__file__))


def _ref_path():
    for d in (os.path.join(os.path.dirname(BASE_DIR), "data"),
              os.path.join(BASE_DIR, "data")):
        p = os.path.join(d, "reference.pkl")
        if os.path.exists(p):
            return p
    return os.path.join(os.path.dirname(BASE_DIR), "data", "reference.pkl")


REF_PKL = _ref_path()
TARGET_PAGES = int(os.environ.get("REPORT_TARGET_PAGES", "100"))


def run_report_job(enc_payload, user_email, display_name="", target_pages=None):
    """Decrypt -> annotate -> render -> email. Returns the variant count.

    `enc_payload` is the encrypted genome (bytes) from crypto.encrypt(); we
    decrypt it here so plaintext exists only in this worker's memory. Safe to
    retry: it has no side effects until the (idempotent enough) email send.
    """
    target_pages = target_pages or TARGET_PAGES
    genome_text = crypto.decrypt(enc_payload).decode("utf-8", errors="replace")
    try:
        bundle = load_reference(REF_PKL)
        records, stats = annotate_genome(genome_text, ref_path=REF_PKL)
        log.info("annotated %d variants", stats["annotated"])

        pdf = build_pdf(records, user_display_name=display_name,
                        target_pages=target_pages, ld_map=bundle.get("ld"))

        if send_report_email and user_email:
            send_report_email(to_address=user_email,
                              user_display_name=display_name,
                              pdf_bytes=pdf,
                              pdf_filename="genodynlabs_report.pdf")
        return stats["annotated"]
    finally:
        # drop references so the GC can reclaim the genome/records promptly
        genome_text = None
