"""
genodyn_report.py  --  Condensed, topic-organized PDF report.

Drop-in alternative to pdf_builder.build_pdf_from_csv(). Same input
(snpedia_data.csv produced by generate_variant_report.py) and same engine
(WeasyPrint), but lays the data out as the reference report does: pastel
topic banners, an "article" sub-header per gene/trait, and dense 5-column
tables (Gene | RS ID | Effect Allele | Your Genotype | Notes) with the
"Your Genotype" cell highlighted.

Content comes from the user's own data:
  * snpedia_data.csv  -> gene, rsid, your genotype, magnitude, repute, summary
  * trait_df.csv      -> effect (risk) allele + named GWAS trait + odds ratio
  * equilibrium_df.csv-> r^2, used to collapse redundant SNPs within an article

We do NOT reproduce any third party's curated notes or article groupings.

Public entry point:
    build_pdf(snpedia_csv, trait_csv, eq_csv, user_display_name="",
              min_magnitude=0.0, max_rows_per_article=14,
              drop_neutral_zero=False, require_note=True) -> bytes
"""
import re
import csv
import math
from io import BytesIO
from collections import defaultdict

import pandas as pd
from weasyprint import HTML, CSS

try:
    from . import topics            # when placed inside the `services` package
except ImportError:
    import topics                   # when run standalone (e.g. CLI/testing)

EMDASH = "\u2014"


# ---------------------------------------------------------------- text utils

def _clean(text, limit=190):
    """Strip SNPedia markup and squeeze a summary down to one tidy line."""
    if text is None:
        return ""
    s = str(text)
    if s.lower() in ("nan", "none"):
        return ""
    s = re.sub(r"\[\[([^\]]+)\]\]", r"\1", s)        # [[Link]] -> Link
    s = s.replace("$", " ")                            # SNPedia newline marker
    s = re.sub(r"\{\{[^}]*\}\}", " ", s)               # leftover templates
    s = re.sub(r"\s+", " ", s).strip(" .;,-")
    if len(s) > limit:
        cut = s[:limit].rsplit(" ", 1)[0]
        s = cut + "\u2026"
    return s


def _esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            if s is not None else "")


def _norm_geno(g):
    g = (str(g) if g is not None else "").strip().upper()
    if g in ("", "NAN", "NONE", "--", "..", "00"):
        return "--"
    return g


def _to_float(x):
    try:
        v = float(x)
        return v if not math.isnan(v) else None
    except (TypeError, ValueError):
        return None


# ------------------------------------------------------- trait_df enrichment

def _build_trait_index(trait_csv):
    """rsid -> {effect_allele, trait, or_txt, n_traits} for the strongest hit."""
    try:
        td = pd.read_csv(trait_csv)
    except Exception:
        return {}
    idx = {}
    for rsid, grp in td.groupby("rsid"):
        # strongest association = smallest p-value
        grp = grp.copy()
        grp["_p"] = grp["pval"].apply(_to_float)
        grp = grp.sort_values("_p", na_position="last")
        best = grp.iloc[0]
        allele = str(best.get("risk_allele", "") or "").strip().upper()
        if allele in ("", "NR", "NAN", "NONE"):
            allele = ""
        orv = _to_float(best.get("OR"))
        rng = str(best.get("range", "") or "").strip()
        or_txt = ""
        if orv:
            or_txt = f"OR {orv:.2f}"
            if rng and rng.upper() not in ("NR", "NAN", ""):
                or_txt += f" {rng}"
        idx[str(rsid)] = {
            "effect_allele": allele,
            "trait": str(best.get("trait", "") or "").strip(),
            "or_txt": or_txt,
            "n_traits": grp["trait"].nunique(),
        }
    return idx


# -------------------------------------------------------------- LD dedup

def _build_ld_index(eq_csv, threshold=0.8):
    """rsid -> set(rsids tightly correlated with it)."""
    ld = defaultdict(set)
    try:
        with open(eq_csv, newline="", encoding="utf-8", errors="replace") as f:
            for row in csv.DictReader(f):
                r2 = _to_float(row.get("r2"))
                if r2 is not None and r2 >= threshold:
                    a, b = row.get("rsid1"), row.get("rsid2")
                    if a and b:
                        ld[a].add(b)
                        ld[b].add(a)
    except Exception:
        pass
    return ld


def _dedup_article(rows, ld):
    """Within one article keep the highest-magnitude SNP of each LD cluster.

    A row that carries a note is never collapsed away — LD-dedup only removes
    redundant rows that have nothing to say. This keeps the dedup from hiding
    informative variants (e.g. two MTHFD1 SNPs with distinct GWAS notes) just
    because they happen to be in linkage disequilibrium.
    """
    rows = sorted(rows, key=lambda r: r["mag"], reverse=True)
    kept, kept_ids = [], set()
    for r in rows:
        has_note = bool((r.get("note") or "").strip())
        if not has_note and any(o in kept_ids for o in ld.get(r["rsid"], ())):
            continue
        kept.append(r)
        kept_ids.add(r["rsid"])
    return kept


# -------------------------------------------------------- genotype highlight

def _highlight(effect_allele, genotype, repute):
    """CSS class for the 'Your Genotype' cell, mirroring the reference look:
    highlight only genotypes worth noticing (you carry the effect allele, or
    SNPedia flags the genotype as notable). Common/protective genotypes are
    left unhighlighted to keep the page readable."""
    g = genotype
    if g == "--":
        return ""
    if effect_allele:
        dose = g.count(effect_allele)
        if dose >= 2:
            return "hl-orange"   # homozygous for effect allele
        if dose == 1:
            return "hl-yellow"   # one copy
        return ""
    if (repute or "").strip().lower() == "bad":
        return "hl-orange"       # notable genotype, no effect allele on file
    return ""


# ----------------------------------------------------------- assemble rows

def _load_rows(snpedia_csv, trait_idx, min_magnitude, drop_neutral_zero,
               require_note):
    rows = []
    with open(snpedia_csv, newline="", encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f):
            rsid = (r.get("rsid") or "").strip()
            if not rsid:
                continue
            mag = _to_float(r.get("mag")) or 0.0
            repute = (r.get("repute") or "neutral").strip() or "neutral"
            gene = (r.get("gene") or "").strip()
            geno = _norm_geno(r.get("geno"))
            summary = _clean(r.get("summary"))
            # SNPedia array-QC notes carry no biological meaning -> drop them
            if re.match(r"(?i)^(common (in|on)\b|normal$|common/normal$|common$)", summary):
                summary = ""

            tinfo = trait_idx.get(rsid, {})
            effect = tinfo.get("effect_allele", "")
            trait = tinfo.get("trait", "")
            or_txt = tinfo.get("or_txt", "")

            # Note: prefer the SNPedia summary; fold in the GWAS trait/OR.
            note = summary
            if trait:
                gwas = trait if not or_txt else f"{trait} ({or_txt})"
                note = f"{gwas}. {note}" if note else gwas
            if not note:
                note = _clean(r.get("geno_description")) or _clean(r.get("description"))

            hl = _highlight(effect, geno, repute)

            topic = topics.classify(gene=gene, trait=trait)
            article = (trait.title() if trait else (gene or "Other variants"))
            article_kind = "trait" if trait else "gene"

            # Quality filter: a row earns its place only if it actually says
            # something — it has a note, a named trait, an effect allele, or a
            # flagged genotype. This is what drops the "Good"/blank filler rows
            # that drop_neutral_zero misses (repute Good but no note), taking
            # the report from ~250 pages down to ~90.
            if require_note and not (note or trait or effect or hl):
                continue

            is_neutral_zero = (mag == 0.0 and repute.lower() == "neutral" and not trait)
            if drop_neutral_zero and is_neutral_zero:
                continue
            if mag < min_magnitude and not trait:
                continue

            rows.append({
                "rsid": rsid, "gene": gene or "\u2014", "effect": effect,
                "geno": geno, "note": note, "mag": mag, "repute": repute,
                "topic": topic, "article": article, "article_kind": article_kind,
                "hl": hl,
            })
    return rows


# --------------------------------------------------------------- rendering

CSS_TEXT = """
@page {
    size: letter; margin: 0.6in 0.6in 0.75in 0.6in;
    @bottom-center { content: "Page " counter(page) " of " counter(pages);
                     font-size: 8pt; color: #999; }
    @bottom-left  { content: "GenoDynLabs \u2014 Genomic Report";
                     font-size: 7.5pt; color: #bbb; }
}
* { box-sizing: border-box; }
body { font-family: "DejaVu Sans", "Liberation Sans", Helvetica, Arial, sans-serif;
       font-size: 9pt; color: #1f1f1f; line-height: 1.35; }

/* cover */
.cover { text-align: center; padding-top: 2.2in; page-break-after: always; }
.cover .brand { font-size: 13pt; letter-spacing: 4pt; color: #6366f1; font-weight: 700; }
.cover h1 { font-size: 30pt; margin: 0.35in 0 0.12in; color: #111; }
.cover .sub { font-size: 12pt; color: #666; }
.cover .meta { font-size: 10pt; color: #999; margin-top: 0.8in; }

/* summary */
.summary { background: #f5f5fb; border-left: 3pt solid #6366f1;
           padding: 10pt 14pt; margin: 0 0 18pt; font-size: 9.5pt; }
.summary h2 { margin: 0 0 6pt; font-size: 13pt; color: #4338ca; }
.summary table { width: 100%; border: none; }
.summary td { border: none; padding: 2pt 0; }
.summary .big { font-size: 15pt; font-weight: 700; color: #111; }

/* topic banner */
.topic { margin: 20pt 0 6pt; padding: 8pt 14pt; border-radius: 3pt;
         font-size: 15pt; font-weight: 600; letter-spacing: 0.5pt;
         page-break-after: avoid; }

/* article header */
h3.article { font-size: 10pt; color: #1f7a4d; margin: 12pt 0 3pt;
             padding-bottom: 1pt; border-bottom: 1pt solid #cfe6da;
             page-break-after: avoid; }
h3.article .meta { color: #999; font-weight: 400; font-size: 8.5pt; }

/* data table */
table.snp { width: 100%; border-collapse: collapse; font-size: 8.3pt;
            page-break-inside: auto; }
table.snp thead { display: table-header-group; }
table.snp th { background: #555; color: #fff; text-align: left;
               padding: 3.5pt 6pt; font-weight: 600; font-size: 8.3pt; }
table.snp td { padding: 3pt 6pt; border: 0.5pt solid #e4e4e4; vertical-align: top; }
table.snp tbody tr:nth-child(even) td { background: #f4f4f4; }
table.snp tr { page-break-inside: avoid; }
.c-gene   { width: 11%; }
.c-rsid   { width: 13%; font-family: "DejaVu Sans Mono", monospace; font-size: 7.8pt; }
.c-effect { width: 8%;  text-align: center; }
.c-geno   { width: 11%; text-align: center; font-weight: 700; }
.c-note   { width: 57%; }

.hl-orange { background: #ffcf9e !important; }
.hl-yellow { background: #fff59d !important; }
.hl-bad    { background: #f8d7d2 !important; }
.hl-good   { background: #d9efce !important; }

.disclaimer { margin-top: 26pt; padding: 12pt; background: #fffbea;
              border: 1px solid #f0e2a0; border-radius: 3pt; font-size: 8pt;
              color: #6b5400; }
"""


def _render_html(grouped, totals, user_display_name):
    meta_line = (f"Prepared for {_esc(user_display_name)}"
                 if user_display_name else "")
    p = ['<!doctype html><html><head><meta charset="utf-8"></head><body>']

    # cover
    p += [
        '<div class="cover">',
        '  <div class="brand">GENODYNLABS</div>',
        '  <h1>Genomic Variant Report</h1>',
        '  <div class="sub">Your SNPs, organized by health topic</div>',
        f'  <div class="meta">{meta_line}</div>',
        '</div>',
    ]

    # summary
    p += [
        '<div class="summary"><h2>Summary</h2><table>',
        f'<tr><td class="big">{totals["variants"]:,}</td>'
        f'<td class="big">{totals["with_trait"]:,}</td>'
        f'<td class="big">{totals["flagged"]:,}</td>'
        f'<td class="big">{totals["topics"]}</td></tr>',
        '<tr><td>annotated variants</td><td>with a known disease association</td>'
        '<td>flagged genotypes</td><td>topic areas covered</td></tr>',
        '</table>',
        '<p style="margin:8pt 0 0;">Variants are grouped by physiological system, '
        'then by gene or trait. A highlighted <b>Your Genotype</b> cell means you '
        'carry the effect allele (orange = two copies, yellow = one) or a genotype '
        'SNPedia flags as notable. Data courtesy of SNPedia and the GWAS Catalog.</p>',
        '</div>',
    ]

    # topic sections
    for topic_key in topics.TOPIC_ORDER:
        if topic_key not in grouped:
            continue
        m = topics.TOPIC_META[topic_key]
        p.append(f'<div class="topic" style="background:{m["bg"]};color:{m["fg"]};">'
                 f'{_esc(m["label"])}</div>')
        articles = grouped[topic_key]
        def _akey(a):
            if a.startswith("Additional variants"):
                return (1, 0)               # always last
            return (0, -max(r["mag"] for r in articles[a]))
        for article in sorted(articles, key=_akey):
            rows = articles[article]
            p.append(f'<h3 class="article">{_esc(article)}'
                     f' <span class="meta">({len(rows)} variant'
                     f'{"s" if len(rows) != 1 else ""})</span></h3>')
            p.append('<table class="snp"><thead><tr>'
                     '<th class="c-gene">Gene</th><th class="c-rsid">RS ID</th>'
                     '<th class="c-effect">Effect</th>'
                     '<th class="c-geno">Your Genotype</th>'
                     '<th class="c-note">Notes</th></tr></thead><tbody>')
            for r in sorted(rows, key=lambda x: x["mag"], reverse=True):
                effect_disp = _esc(r["effect"]) if r["effect"] else EMDASH
                p.append(
                    "<tr>"
                    f'<td class="c-gene">{_esc(r["gene"])}</td>'
                    f'<td class="c-rsid">{_esc(r["rsid"])}</td>'
                    f'<td class="c-effect">{effect_disp}</td>'
                    f'<td class="c-geno {r["hl"]}">{_esc(r["geno"])}</td>'
                    f'<td class="c-note">{_esc(r["note"])}</td>'
                    "</tr>"
                )
            p.append('</tbody></table>')

    p.append(
        '<div class="disclaimer"><b>Important:</b> This report is for educational '
        'and informational purposes only. It is not medical advice and is not '
        'intended to diagnose, treat, or prevent any disease. Genotype calls from '
        'consumer arrays can contain errors; clinically important findings should be '
        'confirmed with a validated test. Always consult a qualified healthcare '
        'professional before acting on genetic information.</div>')
    p.append('</body></html>')
    return "\n".join(p)


# --------------------------------------------------------------- entry point

def build_pdf(snpedia_csv, trait_csv, eq_csv, user_display_name="",
              min_magnitude=0.0, max_rows_per_article=14,
              drop_neutral_zero=False, ld_threshold=0.8, require_note=True):
    trait_idx = _build_trait_index(trait_csv)
    ld = _build_ld_index(eq_csv, threshold=ld_threshold)

    rows = _load_rows(snpedia_csv, trait_idx, min_magnitude, drop_neutral_zero,
                      require_note)

    # group: topic -> article -> [rows]
    grouped = defaultdict(lambda: defaultdict(list))
    for r in rows:
        grouped[r["topic"]][r["article"]].append(r)

    # per article: LD dedup + cap; then merge singleton gene-articles
    for topic_key, articles in list(grouped.items()):
        for article, arows in list(articles.items()):
            deduped = _dedup_article(arows, ld)
            # Hard cap: keep the highest-magnitude rows up to the limit. Rows
            # past the cap are dropped even if they have notes, by design, to
            # keep long articles readable.
            articles[article] = sorted(deduped, key=lambda x: x["mag"],
                                       reverse=True)[:max_rows_per_article]
        keep, extras = {}, []
        for article, arows in articles.items():
            if arows[0]["article_kind"] == "trait" or len(arows) >= 3:
                keep[article] = arows
            else:
                extras.extend(arows)
        if extras:
            keep["Additional variants in this category"] = sorted(
                extras, key=lambda x: x["mag"], reverse=True)[:max(40, max_rows_per_article)]
        grouped[topic_key] = keep

    flagged = with_trait = final_count = 0
    for topic_key, articles in grouped.items():
        for article, arows in articles.items():
            for r in arows:
                final_count += 1
                if r["hl"] in ("hl-orange", "hl-bad"):
                    flagged += 1
                if r["effect"]:
                    with_trait += 1

    totals = {
        "variants": final_count,
        "with_trait": with_trait,
        "flagged": flagged,
        "topics": len([k for k in grouped if grouped[k]]),
    }

    html = _render_html(grouped, totals, user_display_name)
    buf = BytesIO()
    HTML(string=html).write_pdf(buf, stylesheets=[CSS(string=CSS_TEXT)])
    return buf.getvalue()


if __name__ == "__main__":
    import sys
    out = build_pdf(sys.argv[1],
                    sys.argv[2] if len(sys.argv) > 2 else "trait_df.csv",
                    sys.argv[3] if len(sys.argv) > 3 else "equilibrium_df.csv",
                    user_display_name="Sample User")
    with open("report.pdf", "wb") as f:
        f.write(out)
    print(f"wrote report.pdf ({len(out):,} bytes)")
