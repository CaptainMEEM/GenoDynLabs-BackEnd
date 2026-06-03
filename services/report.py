"""
report.py  --  Condensed, topic-organized PDF report (rich-data edition).

The LAYOUT is unchanged from the version you liked: a cover, a summary band,
pastel topic banners, an "article" sub-header per gene/trait, and dense
5-column tables (Gene | RS ID | Effect | Your Genotype | Notes) with the
"Your Genotype" cell highlighted.

What changed is the DATA feeding it. Instead of reading the thin snpedia_data.csv
and a mostly-empty `summary`, this consumes the fully-annotated records produced
by annotate.annotate_genome(): every row already has a real note (notes.py), a
function-first article label (gene_labels.py), an effect allele, dosage, and a
highlight class.

Public entry point:
    build_pdf(records, user_display_name="", target_pages=100) -> bytes
"""
from io import BytesIO
from collections import defaultdict

from weasyprint import HTML, CSS

try:
    from . import topics
except ImportError:
    import topics

EMDASH = "\u2014"


def _esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            if s is not None else "")


# ---------------------------------------------------------- article keying

def _article_for(rec):
    """The article (sub-section) a record belongs to. Curated function labels
    win (e.g. 'Vitamin B12 Status (FUT2)'); otherwise a named trait; otherwise
    the gene; otherwise a catch-all."""
    if rec.get("label_known"):
        return rec["label"]                      # curated, function-first
    if rec.get("primary_trait"):
        # Title-case the trait but keep the gene suffix for traceability
        t = rec["primary_trait"].strip()
        t = t[:1].upper() + t[1:]
        g = rec.get("gene") or ""
        return f"{t} ({g})" if g else t
    return rec.get("gene") or "Other variants"


# ------------------------------------------------------------ dedup by LD

def _interest(rec):
    """Ranking score: flagged + has-OR + magnitude + has-real-note. Drives both
    LD-dedup (keep the most informative of a linked cluster) and the page cap."""
    score = float(rec.get("mag") or 0.0)
    if rec.get("hl") in ("hl-orange", "hl-bad"):
        score += 3.0
    if rec.get("hl") == "hl-yellow":
        score += 1.0
    if rec.get("max_or"):
        score += min(abs(rec["max_or"] - 1.0) * 2.0, 3.0)
    if rec.get("label_known") or rec.get("primary_trait"):
        score += 1.0
    return score


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

.cover { text-align: center; padding-top: 2.2in; page-break-after: always; }
.cover .brand { font-size: 13pt; letter-spacing: 4pt; color: #6366f1; font-weight: 700; }
.cover h1 { font-size: 30pt; margin: 0.35in 0 0.12in; color: #111; }
.cover .sub { font-size: 12pt; color: #666; }
.cover .meta { font-size: 10pt; color: #999; margin-top: 0.8in; }

.summary { background: #f5f5fb; border-left: 3pt solid #6366f1;
           padding: 10pt 14pt; margin: 0 0 18pt; font-size: 9.5pt; }
.summary h2 { margin: 0 0 6pt; font-size: 13pt; color: #4338ca; }
.summary table { width: 100%; border: none; }
.summary td { border: none; padding: 2pt 0; }
.summary .big { font-size: 15pt; font-weight: 700; color: #111; }

.topic { margin: 20pt 0 6pt; padding: 8pt 14pt; border-radius: 3pt;
         font-size: 15pt; font-weight: 600; letter-spacing: 0.5pt;
         page-break-after: avoid; }

h3.article { font-size: 10pt; color: #1f7a4d; margin: 12pt 0 3pt;
             padding-bottom: 1pt; border-bottom: 1pt solid #cfe6da;
             page-break-after: avoid; }
h3.article .meta { color: #999; font-weight: 400; font-size: 8.5pt; }

table.snp { width: 100%; border-collapse: collapse; font-size: 8.3pt;
            page-break-inside: auto; }
table.snp thead { display: table-header-group; }
table.snp th { background: #555; color: #fff; text-align: left;
               padding: 3.5pt 6pt; font-weight: 600; font-size: 8.3pt; }
table.snp td { padding: 3pt 6pt; border: 0.5pt solid #e4e4e4; vertical-align: top; }
table.snp tbody tr:nth-child(even) td { background: #f4f4f4; }
table.snp tr { page-break-inside: avoid; }
.c-gene   { width: 11%; }
.c-rsid   { width: 12%; font-family: "DejaVu Sans Mono", monospace; font-size: 7.8pt; }
.c-effect { width: 7%;  text-align: center; }
.c-geno   { width: 10%; text-align: center; font-weight: 700; }
.c-note   { width: 60%; }

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
    p += [
        '<div class="cover">',
        '  <div class="brand">GENODYNLABS</div>',
        '  <h1>Genomic Variant Report</h1>',
        '  <div class="sub">Your SNPs, organized by health topic</div>',
        f'  <div class="meta">{meta_line}</div>',
        '</div>',
    ]
    p += [
        '<div class="summary"><h2>Summary</h2><table>',
        f'<tr><td class="big">{totals["variants"]:,}</td>'
        f'<td class="big">{totals["with_trait"]:,}</td>'
        f'<td class="big">{totals["flagged"]:,}</td>'
        f'<td class="big">{totals["topics"]}</td></tr>',
        '<tr><td>annotated variants</td><td>with a known trait association</td>'
        '<td>flagged genotypes</td><td>topic areas covered</td></tr>',
        '</table>',
        '<p style="margin:8pt 0 0;">Variants are grouped by physiological system, '
        'then by the gene\'s primary function. A highlighted <b>Your Genotype</b> '
        'cell means you carry the effect allele (orange = two copies, yellow = one) '
        'or a genotype SNPedia flags as notable. Notes combine SNPedia genotype '
        'summaries with GWAS Catalog associations (odds ratios shown where known). '
        'Data courtesy of SNPedia and the NHGRI-EBI GWAS Catalog.</p>',
        (f'<p style="margin:6pt 0 0;color:#666;">This report highlights your '
         f'{totals["variants"]:,} most informative variants. '
         f'{totals["truncated"]:,} additional lower-signal matches were '
         f'summarized but not listed individually.</p>'
         if totals["truncated"] > 0 else ""),
        '</div>',
    ]

    for topic_key in topics.TOPIC_ORDER:
        if topic_key not in grouped or not grouped[topic_key]:
            continue
        m = topics.TOPIC_META[topic_key]
        p.append(f'<div class="topic" style="background:{m["bg"]};color:{m["fg"]};">'
                 f'{_esc(m["label"])}</div>')
        articles = grouped[topic_key]

        def _akey(a):
            if a.startswith("Additional variants"):
                return (1, 0)
            return (0, -max(_interest(r) for r in articles[a]))

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
            for r in sorted(rows, key=_interest, reverse=True):
                effect_disp = _esc(r["effect_allele"]) if r["effect_allele"] else EMDASH
                p.append(
                    "<tr>"
                    f'<td class="c-gene">{_esc(r.get("gene") or EMDASH)}</td>'
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
        'intended to diagnose, treat, or prevent any disease. Many traits are '
        'polygenic, and the effect of any single SNP is usually small; per-SNP '
        'summaries from SNPedia can be incomplete or out of date. Genotype calls '
        'from consumer arrays can contain errors; clinically important findings '
        'should be confirmed with a validated clinical test. Always consult a '
        'qualified healthcare professional before acting on genetic information.'
        '</div>')
    p.append('</body></html>')
    return "\n".join(p)


# --------------------------------------------------------------- entry point

def build_pdf(records, user_display_name="", target_pages=100,
              max_rows_per_article=18, ld_map=None, rows_per_page=12):
    """Group, de-duplicate, cap, and render the annotated records to PDF bytes.

    A GLOBAL row budget (target_pages * rows_per_page) keeps the report near the
    requested length no matter how dense the input is: a sequenced genome that
    matches 12k SNPedia entries and a sparse consumer chip that matches 4k both
    land near `target_pages`. We keep the highest-`_interest` variants first, so
    flagged genotypes and strong trait associations are never the ones dropped.
    Every row that survives still carries a full note.

    `records` come from annotate.annotate_genome(). `ld_map` (rsid -> [rsids])
    is optional; pass bundle["ld"] to collapse linked SNPs within an article.
    """
    ld_map = ld_map or {}

    total_found = len(records)
    budget = max(400, int(target_pages * rows_per_page))

    # group ALL records: topic -> article -> [records]
    grouped_all = defaultdict(lambda: defaultdict(list))
    for r in records:
        topic = topics.classify(gene=r.get("gene"), trait=r.get("primary_trait"))
        grouped_all[topic][_article_for(r)].append(r)

    # de-dup + cap each article, then collapse the weak single-row tail into a
    # per-topic "Additional" table so the report reads as multi-row sections.
    grouped = defaultdict(dict)
    for topic_key, articles in grouped_all.items():
        capped = {}
        for article, arows in articles.items():
            kept, kept_ids = [], set()
            for r in sorted(arows, key=_interest, reverse=True):
                links = ld_map.get(r["rsid"], ())
                informative = (r.get("hl") in ("hl-orange", "hl-yellow")
                               or r.get("max_or") or r.get("label_known"))
                if not informative and any(o in kept_ids for o in links):
                    continue
                kept.append(r)
                kept_ids.add(r["rsid"])
            capped[article] = kept[:max_rows_per_article]

        keep, extras = {}, []
        for article, arows in capped.items():
            standout = any(r.get("hl") == "hl-orange"
                           or (r.get("mag") or 0) >= 2.5 for r in arows)
            # own header if: 2+ variants, OR a genuinely notable single finding
            if len(arows) >= 2 or standout:
                keep[article] = arows
            else:
                extras.extend(arows)
        if extras:
            keep["Additional variants in this category"] = sorted(
                extras, key=_interest, reverse=True)[:120]
        grouped[topic_key] = keep

    # ---- article-first global budget. Phase 1: select the strongest NAMED
    # articles (all their rows together) until the row budget fills. Phase 2:
    # spend any leftover budget on the per-topic "Additional" tables, trimmed
    # to fit. This keeps page count near target with dense, multi-row sections.
    named, additional = [], []
    for topic_key, articles in grouped.items():
        for article, arows in articles.items():
            if article.startswith("Additional"):
                additional.append((topic_key, article, arows))
            else:
                named.append((topic_key, article, arows,
                              max(_interest(r) for r in arows)))
    named.sort(key=lambda x: x[3], reverse=True)

    chosen = defaultdict(dict)
    rows_used = 0
    for topic_key, article, arows, _score in named:
        if rows_used >= budget:
            break
        chosen[topic_key][article] = arows
        rows_used += len(arows)

    # Phase 2: fill remaining room with Additional tables (strongest topics
    # first), trimming the last one so we don't blow the budget.
    additional.sort(key=lambda x: max(_interest(r) for r in x[2]), reverse=True)
    for topic_key, article, arows in additional:
        room = budget - rows_used
        if room <= 0:
            break
        take = arows[:room]
        chosen[topic_key][article] = take
        rows_used += len(take)

    grouped = chosen
    truncated = total_found - rows_used

    # totals
    flagged = with_trait = final_count = 0
    for articles in grouped.values():
        for arows in articles.values():
            for r in arows:
                final_count += 1
                if r.get("hl") in ("hl-orange", "hl-bad"):
                    flagged += 1
                if r.get("primary_trait") or r.get("effect_allele"):
                    with_trait += 1

    totals = {
        "variants": final_count,
        "with_trait": with_trait,
        "flagged": flagged,
        "topics": len([k for k in grouped if grouped[k]]),
        "total_found": total_found,
        "truncated": truncated,
    }

    html = _render_html(grouped, totals, user_display_name)
    buf = BytesIO()
    HTML(string=html).write_pdf(buf, stylesheets=[CSS(string=CSS_TEXT)])
    return buf.getvalue()
