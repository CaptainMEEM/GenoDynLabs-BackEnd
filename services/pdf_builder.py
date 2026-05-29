"""
Build a clean PDF from the annotated CSV.

We don't render the interactive variant_report.html directly because it relies
on JavaScript at runtime (the data is base64-encoded into a JS variable, then
DataTables renders the rows in the browser). WeasyPrint doesn't execute JS,
so we build a static PDF directly from the CSV instead — same data, just
laid out for paper.
"""
import csv
from io import BytesIO
from weasyprint import HTML, CSS


PDF_CSS = """
@page {
    size: letter;
    margin: 0.75in;
    @bottom-center {
        content: "Page " counter(page) " of " counter(pages);
        font-size: 9pt;
        color: #666;
    }
}
body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
    font-size: 10pt;
    color: #1a1a1a;
    line-height: 1.4;
}
.cover {
    text-align: center;
    padding-top: 2in;
    page-break-after: always;
}
.cover .brand {
    font-size: 14pt;
    letter-spacing: 3pt;
    color: #6366f1;
    font-weight: 700;
}
.cover h1 {
    font-size: 28pt;
    margin: 0.4in 0 0.15in;
    color: #111;
}
.cover .subtitle {
    font-size: 12pt;
    color: #555;
    margin-bottom: 0.8in;
}
.cover .meta {
    font-size: 10pt;
    color: #888;
}
h2 {
    font-size: 14pt;
    color: #4338ca;
    border-bottom: 2px solid #6366f1;
    padding-bottom: 4pt;
    margin-top: 24pt;
}
.summary-box {
    background: #f5f5f9;
    border-left: 3pt solid #6366f1;
    padding: 10pt 14pt;
    margin: 12pt 0 20pt;
    font-size: 10pt;
}
.variant {
    border: 1px solid #e0e0e8;
    border-radius: 4pt;
    padding: 8pt 10pt;
    margin-bottom: 8pt;
    page-break-inside: avoid;
}
.variant-header {
    display: flex;
    justify-content: space-between;
    font-weight: 600;
    margin-bottom: 4pt;
}
.rsid { color: #4338ca; font-family: "SF Mono", Menlo, Consolas, monospace; }
.gene { color: #555; font-size: 9pt; }
.repute-good  { color: #059669; font-weight: 600; }
.repute-bad   { color: #dc2626; font-weight: 600; }
.repute-neutral { color: #6b7280; font-weight: 600; }
.summary { font-size: 10pt; margin: 4pt 0; }
.geno-line { font-size: 9pt; color: #555; }
.disclaimer {
    margin-top: 40pt;
    padding: 14pt;
    background: #fffbea;
    border: 1px solid #f5e89c;
    border-radius: 4pt;
    font-size: 9pt;
    color: #6b5400;
}
"""


def _repute_class(repute: str) -> str:
    r = (repute or "").strip().lower()
    if r == "good":
        return "repute-good"
    if r == "bad":
        return "repute-bad"
    return "repute-neutral"


def _esc(s) -> str:
    if s is None:
        return ""
    return (str(s)
            .replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;"))


def build_pdf_from_csv(csv_path: str, user_display_name: str = "") -> bytes:
    """
    Read snpedia_data.csv and produce a PDF report. Returns the PDF bytes.
    """
    rows = []
    with open(csv_path, newline="", encoding="utf-8", errors="replace") as f:
        for row in csv.DictReader(f):
            # parse magnitude as float for sorting (blank = 0)
            try:
                mag = float(row.get("mag", "") or 0)
            except ValueError:
                mag = 0.0
            row["_mag"] = mag
            rows.append(row)

    # Sort by magnitude descending — most significant first
    rows.sort(key=lambda r: r["_mag"], reverse=True)
    total = len(rows)
    significant = [r for r in rows if r["_mag"] > 0]

    # Build the HTML
    name_line = f"Prepared for {_esc(user_display_name)}" if user_display_name else ""

    parts = [
        "<!doctype html><html><head><meta charset='utf-8'></head><body>",
        # Cover page
        "<div class='cover'>",
        "  <div class='brand'>GENODYNLABS</div>",
        "  <h1>Genomic Variant Report</h1>",
        "  <div class='subtitle'>Annotated SNP analysis based on SNPedia</div>",
        f"  <div class='meta'>{name_line}</div>",
        "</div>",
        # Summary
        "<h2>Summary</h2>",
        "<div class='summary-box'>",
        f"  <p><strong>Total annotated variants:</strong> {total:,}</p>",
        f"  <p><strong>Variants with reported significance (magnitude &gt; 0):</strong> {len(significant):,}</p>",
        "  <p>The variants below are ordered from highest magnitude (most studied or significant)",
        "  to lowest. Magnitude is SNPedia's measure of how interesting the variant is.</p>",
        "</div>",
        "<h2>Variants</h2>",
    ]

    for r in rows:
        rsid = _esc(r.get("rsid"))
        gene = _esc(r.get("gene"))
        geno = _esc(r.get("geno"))
        repute = _esc(r.get("repute") or "neutral")
        summary = _esc(r.get("summary"))
        mag = r["_mag"]
        chrom = _esc(r.get("chr"))
        pos = _esc(r.get("pos"))

        parts.append("<div class='variant'>")
        parts.append("  <div class='variant-header'>")
        parts.append(f"    <span class='rsid'>{rsid}</span>")
        parts.append(f"    <span class='gene'>Gene: {gene or '—'} · Chr {chrom}:{pos} · Mag {mag:g}</span>")
        parts.append("  </div>")
        parts.append(f"  <div class='geno-line'>Your genotype: <strong>{geno}</strong> · "
                     f"<span class='{_repute_class(repute)}'>{repute.title()}</span></div>")
        if summary:
            parts.append(f"  <div class='summary'>{summary}</div>")
        parts.append("</div>")

    parts.append(
        "<div class='disclaimer'>"
        "<strong>Important:</strong> This report is for educational and informational "
        "purposes only. It is not medical advice and should not be used to diagnose, "
        "treat, or prevent any disease. Always consult a qualified healthcare professional "
        "before making decisions based on genetic information."
        "</div>"
    )
    parts.append("</body></html>")

    html_str = "\n".join(parts)
    buf = BytesIO()
    HTML(string=html_str).write_pdf(buf, stylesheets=[CSS(string=PDF_CSS)])
    return buf.getvalue()
