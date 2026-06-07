"""
precompute_reference.py  --  Build the rich reference bundle ONCE at deploy time.

This is the single most important file for report quality. The OLD version threw
away ~90% of SNPedia's content (it kept only gene/chr/pos/orientation per SNP and
geno/mag/repute/summary per genotype), which is why ~63% of genotype notes came
out blank or boilerplate.

This version mines EVERYTHING the four SNPedia/GWAS dumps contain and folds it
into one compact, rsid-keyed pickle so the worker never touches the 33 MB raw
files or runs a regex at request time:

  snp_df.csv          (SNP pages)        -> gene(s), chr, pos, orientation, GMAF,
                                            and every {{PMID Auto GWAS}} block
                                            (trait, risk allele, OR, study title)
  geno_df.csv         (genotype pages)   -> per-genotype mag, repute, summary, desc
  trait_df.csv        (GWAS Catalog)     -> risk allele, OR, named trait, population
                                            allele frequencies (per-rsid)
  equilibrium_df.csv  (LD r^2)           -> linkage map for de-duplication

Output schema (reference.pkl):

    {
      "version": "2.0",
      "snps": {
        "rs602662": {
          "rsid": "rs602662",
          "gene": "FUT2",
          "genes": ["FUT2"],
          "chr": "19", "pos": "48703417",
          "orientation": "minus",
          "gmaf": 0.43,
          "gwas": [
            {"trait": "Folate pathway vitamins", "risk": "G",
             "or": 1.20, "title": "...", "pval": "3E-12"}
          ],
          "catalog": [            # from trait_df (GWAS Catalog)
            {"trait": "Vitamin B12 levels", "risk": "G", "or": 1.31,
             "range": "[1.2-1.4]", "pval": "...", "raf": {"European": 0.49, ...}}
          ],
          "options": {            # keyed by SORTED genotype, e.g. "AG"
            "GG": {"mag": 2.0, "repute": "Good",
                   "summary": "Higher vitamin B12 levels", "desc": "..."},
            ...
          }
        },
        ...
      },
      "ld": { "rs1234": ["rs5678", ...], ... }   # r^2 >= LD_THRESHOLD
    }

Run:
    python -m services.precompute_reference \
        data/snp_df.csv data/geno_df.csv data/trait_df.csv \
        data/equilibrium_df.csv data/reference.pkl
"""
import re
import os
import sys
import pickle
from collections import defaultdict

import pandas as pd

LD_THRESHOLD = 0.8          # only keep tight linkage in the pickle (dedup use)
MAX_GWAS_PER_SNP = 6        # cap associations so notes stay readable
MAX_SUMMARY_LEN = 320       # genotype summaries are short; keep them whole-ish


# ---------------------------------------------------------------- field utils
# SNPedia dumps use `$` for newlines and `|` to separate fields inside {{ }}.

def _field(name, text):
    """First value of `name=...` inside a SNPedia template blob."""
    m = re.search(r"%s=([^$|}]+)" % re.escape(name), text or "")
    return m.group(1).strip() if m else ""


def _to_float(x):
    try:
        v = float(x)
        return v if v == v else None        # filter NaN
    except (TypeError, ValueError):
        return None


def _clean_text(text, limit=MAX_SUMMARY_LEN):
    """Strip SNPedia/wiki markup down to a clean human sentence."""
    if not text:
        return ""
    s = re.sub(r"\{\{[^}]*\}\}", " ", text)       # drop template blocks
    s = re.sub(r"\[\[([^\]|]+\|)?([^\]]+)\]\]", r"\2", s)  # [[a|b]] -> b
    s = s.replace("$", " ")
    s = re.sub(r"\s+", " ", s).strip(" .;,-|")
    if len(s) > limit:
        s = s[:limit].rsplit(" ", 1)[0] + "\u2026"
    return s


# ------------------------------------------------------------ snp_df parsing

def _parse_snp_page(text):
    """Return (rsid, record) mined from one SNP wiki page, or (None, None)."""
    rsnum = _field("rsid", text)
    if not rsnum:
        return None, None
    rsid = "rs" + rsnum

    genes_raw = _field("Gene_s", text) or _field("Gene", text)
    genes = [g.strip() for g in re.split(r"[,;]", genes_raw) if g.strip()]
    primary = _field("Gene", text) or (genes[0] if genes else "")

    # Every GWAS association attached to this SNP page.
    gwas = []
    for block in re.findall(r"\{\{PMID Auto GWAS(.*?)\}\}", text, re.DOTALL):
        trait = _field("Trait", block)
        if not trait or trait.lower() in ("none", "nr", ""):
            continue
        risk = _field("RiskAllele", block).upper()
        if risk in ("NR", "NONE", "?"):
            risk = ""
        gwas.append({
            "trait": trait,
            "risk": risk,
            "or": _to_float(_field("OR", block)),
            "title": _field("Title", block),
            "pval": _field("Pval", block),
        })
    # de-dup identical (trait, risk) keeping the most extreme OR; cap length
    seen = {}
    for g in gwas:
        key = (g["trait"].lower(), g["risk"])
        prev = seen.get(key)
        if prev is None or abs((g["or"] or 1) - 1) > abs((prev["or"] or 1) - 1):
            seen[key] = g
    gwas = sorted(seen.values(),
                  key=lambda g: abs((g["or"] or 1) - 1), reverse=True)[:MAX_GWAS_PER_SNP]

    rec = {
        "rsid": rsid,
        "gene": primary,
        "genes": genes,
        "chr": _field("Chromosome", text),
        "pos": _field("position", text),
        "orientation": (_field("StabilizedOrientation", text)
                        or _field("Orientation", text)).lower(),
        "gmaf": _to_float(_field("GMAF", text)),
        "gwas": gwas,
        "catalog": [],          # filled from trait_df later
        "options": {},          # filled from geno_df later
    }
    return rsid, rec


def parse_snps(snp_csv):
    snps = {}
    df = pd.read_csv(snp_csv)
    for text in df["text"].dropna():
        rsid, rec = _parse_snp_page(str(text))
        if rsid:
            snps[rsid] = rec
    return snps


# ----------------------------------------------------------- geno_df parsing

def _parse_genotype_cell(text):
    """Return (rsid, sorted_geno, option) for one genotype page cell."""
    if not text or text == "none":
        return None, None, None
    rsnum = _field("rsid", text)
    a1, a2 = _field("allele1", text), _field("allele2", text)
    geno = "".join(sorted([a1, a2])).lstrip("=").strip()
    if not rsnum or not geno:
        return None, None, None
    summary = _field("summary", text)
    opt = {
        "mag": _to_float(_field("magnitude", text)) or 0.0,
        "repute": (_field("repute", text) or "neutral"),
        "summary": _clean_text(summary) if summary else "",
        # full-page text minus the template noise; a richer fallback note
        "desc": _clean_text(text, limit=MAX_SUMMARY_LEN) if not summary else "",
    }
    return "rs" + rsnum, geno, opt


def attach_genotypes(snps, geno_csv):
    """Fold every genotype option into the matching SNP record.

    A genotype page may reference an rsid that has no SNP page; in that case we
    create a minimal SNP record so the variant is never silently dropped.
    """
    df = pd.read_csv(geno_csv)
    for row in df.itertuples(index=False):
        for cell in (row.gt1, row.gt2, row.gt3):
            rsid, geno, opt = _parse_genotype_cell(str(cell))
            if not rsid:
                continue
            rec = snps.get(rsid)
            if rec is None:
                rec = {"rsid": rsid, "gene": "", "genes": [], "chr": "",
                       "pos": "", "orientation": "", "gmaf": None,
                       "gwas": [], "catalog": [], "options": {}}
                snps[rsid] = rec
            rec["options"][geno] = opt
    return snps


# ---------------------------------------------------------- trait_df parsing

_POP_COLS = ["Global", "European", "African", "AfricanOthers", "AfricanAmerican",
             "Asian", "EastAsian", "OtherAsian", "LatinAmerican1",
             "LatinAmerican2", "SouthAsian", "Other"]


def attach_catalog(snps, trait_csv):
    """Fold GWAS-Catalog rows (trait_df) into each SNP's `catalog` list.

    These carry the most trustworthy, study-backed effect sizes and the
    population allele frequencies the trait report uses. We keep the strongest
    few per rsid (smallest p-value)."""
    td = pd.read_csv(trait_csv)
    td["_p"] = pd.to_numeric(td.get("pval"), errors="coerce")
    for rsid, grp in td.groupby("rsid"):
        rec = snps.get(rsid)
        grp = grp.sort_values("_p", na_position="last")
        cats = []
        seen = set()
        for r in grp.itertuples(index=False):
            trait = str(getattr(r, "trait", "") or "").strip()
            if not trait or trait.lower() in ("nan", "none", "nr"):
                continue
            key = trait.lower()
            if key in seen:
                continue
            seen.add(key)
            raf = {}
            for c in _POP_COLS:
                v = _to_float(getattr(r, c, None))
                if v is not None and 0 < v < 1:
                    raf[c] = round(v, 4)
            risk = str(getattr(r, "risk_allele", "") or "").strip().upper()
            if risk in ("NR", "NONE", "?", "NAN"):
                risk = ""
            cats.append({
                "trait": trait,
                "risk": risk,
                "or": _to_float(getattr(r, "OR", None)),
                "range": str(getattr(r, "range", "") or "").strip(),
                "pval": str(getattr(r, "pval", "") or "").strip(),
                "gene": str(getattr(r, "gene", "") or "").strip(),
                "raf": raf,
            })
            if len(cats) >= 6:
                break
        if rec is None:
            # rsid known to GWAS Catalog but absent from SNPedia: keep a stub so
            # the trait association can still surface if the user is genotyped.
            rec = {"rsid": rsid, "gene": cats[0]["gene"] if cats else "",
                   "genes": [], "chr": "", "pos": "", "orientation": "",
                   "gmaf": None, "gwas": [], "catalog": [], "options": {}}
            snps[rsid] = rec
        rec["catalog"] = cats
        if not rec["gene"] and cats and cats[0]["gene"]:
            rec["gene"] = cats[0]["gene"]
    return snps


# ------------------------------------------------- deana nutrient enrichment
# We pull ONLY nutrient level/status associations from the deana evidence pack
# (a few thousand records) so the bundle stays small and the rest of the report
# is unchanged. Each becomes a catalog-style association on its SNP, tagged with
# the nutrient so the report routes it into the Vitamins & Minerals panel.

# Trait-name -> nutrient (panel label). Matched on the record TITLE, and the
# title must also contain a level/status word, so we only ingest genuine
# circulating-level associations (not unrelated clinical variants).
import re as _re
_NUTRIENT_TRAIT = [
    ("Vitamin A", _re.compile(r"\bretinol\b|beta-?carotene|\bvitamin a\b|carotenoid", _re.I)),
    ("Vitamin B2 (Riboflavin)", _re.compile(r"riboflavin", _re.I)),
    ("Vitamin B6", _re.compile(r"pyridoxal|pyridoxine|\bvitamin b6\b|\bPLP\b", _re.I)),
    ("Vitamin B9 (Folate)", _re.compile(r"\bfolate\b|folic acid", _re.I)),
    ("Vitamin B12", _re.compile(r"\bvitamin b-?12\b|\bb12\b|cobalamin", _re.I)),
    ("Vitamin C", _re.compile(r"ascorb|\bvitamin c\b", _re.I)),
    ("Vitamin D", _re.compile(r"\bvitamin d\b|25-?hydroxyvitamin|calcidiol", _re.I)),
    ("Vitamin E", _re.compile(r"\bvitamin e\b|tocopherol|tocotrienol", _re.I)),
    ("Vitamin K", _re.compile(r"phylloquinone|\bvitamin k\b|menaquinone", _re.I)),
    ("Calcium", _re.compile(r"\bcalcium\b", _re.I)),
    ("Magnesium", _re.compile(r"magnesium", _re.I)),
    ("Sodium", _re.compile(r"\bsodium\b", _re.I)),
    ("Potassium", _re.compile(r"potassium", _re.I)),
    ("Iron", _re.compile(r"\biron\b|ferritin|transferrin saturation|iron status", _re.I)),
    ("Zinc", _re.compile(r"\bzinc\b", _re.I)),
    ("Copper", _re.compile(r"\bcopper\b|ceruloplasmin", _re.I)),
    ("Manganese", _re.compile(r"manganese", _re.I)),
    ("Iodine", _re.compile(r"\biodine\b", _re.I)),
    ("Selenium", _re.compile(r"selenium", _re.I)),
    ("Molybdenum", _re.compile(r"molybden", _re.I)),
]
_LEVEL_WORD = _re.compile(
    r"level|status|concentration|plasma|serum|circulating|biomarker", _re.I)
# avoid a couple of well-known false positives that contain a nutrient word
_TRAIT_EXCLUDE = _re.compile(
    r"coronary artery calc|arterial calc|valve calc|channel blocker", _re.I)


def _nutrient_for_trait(title):
    if not title or _TRAIT_EXCLUDE.search(title) or not _LEVEL_WORD.search(title):
        return None
    for label, rx in _NUTRIENT_TRAIT:
        if rx.search(title):
            return label
    return None


def attach_deana_nutrients(snps, deana_shards_dir, max_assoc_per_snp=6):
    """Fold nutrient level/status associations from the deana evidence pack into
    `snps`. Streams shard JSON files one at a time (build-time only) and keeps
    only nutrient-relevant records, so memory and bundle size stay small.

    Each match: ensures the SNP record exists, tags it with `nutrient` (panel
    routing), and appends a catalog-style association the note engine already
    understands. New nutrient rsids become matchable; everything else is ignored
    so the rest of the report is unchanged."""
    import glob
    import json
    shard_files = sorted(glob.glob(os.path.join(deana_shards_dir, "*.json")))
    if not shard_files:
        print(f"  (no deana shards found in {deana_shards_dir}; skipping)")
        return snps

    added_records = 0
    enriched_snps = set()
    seen = set()                       # (rsid, trait) de-dup
    for sf in shard_files:
        try:
            with open(sf) as fh:
                recs = json.load(fh)
        except Exception:
            continue
        for r in recs:
            title = r.get("title", "")
            nut = _nutrient_for_trait(title)
            if not nut:
                continue
            trait = title.split(" association near")[0].strip()
            risk = (r.get("riskAllele") or "").strip().upper()
            if risk in ("NR", "NONE", "?", "NAN"):
                risk = ""
            genes_raw = r.get("genes") or []
            gene = ""
            for g in genes_raw:
                gene = str(g).replace(" - ", "/").split("/")[0].strip()
                if gene:
                    break
            for rsid in r.get("markerIds", []):
                if not rsid or not rsid.startswith("rs"):
                    continue
                key = (rsid, trait.lower())
                if key in seen:
                    continue
                seen.add(key)
                # Only ADD genuinely new nutrient-level SNPs. If the rsid already
                # exists in our reference it has an established meaning/topic
                # (e.g. APOE = Alzheimer); we never re-tag it, so the rest of the
                # report is unchanged and no variant gets hijacked into the panel.
                if rsid in snps:
                    continue
                rec = {"rsid": rsid, "gene": gene, "genes": genes_raw,
                       "chr": "", "pos": "", "orientation": "", "gmaf": None,
                       "gwas": [], "catalog": [], "options": {},
                       "nutrient": nut}
                snps[rsid] = rec
                added_records += 1
                rec["catalog"].append({
                    "trait": trait, "risk": risk, "or": None,
                    "range": "", "pval": "", "gene": gene, "raf": {},
                    "source": "deana", "evidence": r.get("evidenceLevel", ""),
                    "repute": r.get("repute", ""),
                })
                enriched_snps.add(rsid)
    print(f"  deana nutrient associations: {len(enriched_snps):,} SNPs "
          f"enriched ({added_records:,} newly added)")
    return snps


def attach_deana_general(snps, deana_shards_dir, levels=("high", "moderate")):
    """Fold the BROADER deana evidence pack (clinical + GWAS) into the bundle so
    that unique user uploads get matched and annotated. To stay lean we keep
    only `levels` evidence (default high+moderate, dropping ~80k weak
    'supplementary' records) and only ADD new rsids (never touch existing ones,
    so the curated sections and SNPedia notes are unchanged). Each new SNP gets
    one compact catalog entry carrying the trait and, crucially, the ClinVar
    clinical significance so pathogenic findings can be flagged."""
    import glob
    import json
    levels = set(levels)
    shard_files = sorted(glob.glob(os.path.join(deana_shards_dir, "*.json")))
    if not shard_files:
        print(f"  (no deana shards in {deana_shards_dir}; skipping general merge)")
        return snps
    added = 0
    for sf in shard_files:
        try:
            recs = json.load(open(sf))
        except Exception:
            continue
        for r in recs:
            if r.get("evidenceLevel") not in levels:
                continue
            ids = [m for m in r.get("markerIds", []) if m and m.startswith("rs")]
            if not ids:
                continue
            title = (r.get("title") or "").strip()
            trait = title.split(" association near")[0].strip()[:80]
            risk = (r.get("riskAllele") or "").strip().upper()
            if risk in ("NR", "NONE", "?", "NAN"):
                risk = ""
            gene = ""
            for g in (r.get("genes") or []):
                gene = str(g).replace(" - ", "/").split("/")[0].strip()
                if gene:
                    break
            sig = r.get("clinicalSignificance") or ""
            for rsid in ids:
                if rsid in snps:                 # never overwrite existing
                    continue
                snps[rsid] = {
                    "rsid": rsid, "gene": gene, "genes": ([gene] if gene else []),
                    "chr": "", "pos": "", "orientation": "", "gmaf": None,
                    "gwas": [], "options": {},
                    "catalog": [{"trait": trait, "risk": risk, "or": None,
                                 "range": "", "pval": "", "gene": gene, "raf": {},
                                 "source": "deana", "sig": sig,
                                 "evidence": r.get("evidenceLevel", "")}],
                }
                added += 1
    print(f"  deana general merge ({'+'.join(sorted(levels))}): +{added:,} new SNPs")
    return snps


def ensure_curated(snps, curated_map, baselines):
    """Guarantee every user-specified variant ID is present and matchable. IDs
    already in the bundle keep their data; missing ones get a stub carrying the
    category baseline note (so there is always *some* notation)."""
    added = 0
    for rsid, (label, _topic) in curated_map.items():
        if rsid in snps:
            continue
        note = baselines.get(label, "A user-specified variant.")
        snps[rsid] = {
            "rsid": rsid, "gene": "", "genes": [],
            "chr": "", "pos": "", "orientation": "", "gmaf": None,
            "gwas": [], "options": {},
            "catalog": [{"trait": note, "risk": "", "or": None, "range": "",
                         "pval": "", "gene": "", "raf": {}, "source": "curated"}],
        }
        added += 1
    print(f"  curated variants ensured present: +{added:,} stubs "
          f"({len(curated_map)} total IDs)")
    return snps


# ------------------------------------------------------- equilibrium parsing

def build_ld(eq_csv, threshold=LD_THRESHOLD):
    ld = defaultdict(set)
    df = pd.read_csv(eq_csv)
    for r in df.itertuples(index=False):
        r2 = _to_float(getattr(r, "r2", None))
        if r2 is not None and r2 >= threshold:
            a, b = getattr(r, "rsid1", None), getattr(r, "rsid2", None)
            if a and b:
                ld[a].add(b)
                ld[b].add(a)
    return {k: sorted(v) for k, v in ld.items()}


# --------------------------------------------------------------------- build

def build(snp_csv, geno_csv, trait_csv, eq_csv, out_pkl, deana_dir=None,
          deana_levels=("high", "moderate")):
    print("parsing SNP pages ...")
    snps = parse_snps(snp_csv)
    print(f"  {len(snps):,} SNP pages")

    print("attaching genotype options ...")
    attach_genotypes(snps, geno_csv)

    print("attaching GWAS-Catalog associations ...")
    attach_catalog(snps, trait_csv)

    if deana_dir:
        print("attaching deana nutrient-level associations ...")
        attach_deana_nutrients(snps, deana_dir)
        if deana_levels:
            print("attaching deana general evidence (broad coverage) ...")
            attach_deana_general(snps, deana_dir, levels=deana_levels)

    # Guarantee the user-specified variant IDs are present and routed.
    try:
        from . import curated as _curated
    except ImportError:
        import curated as _curated
    print("ensuring curated variant IDs ...")
    ensure_curated(snps, _curated.CURATED_VARIANTS, _curated.CURATED_BASELINE)

    print("building linkage map ...")
    ld = build_ld(eq_csv)
    print(f"  {len(ld):,} SNPs with tight LD")

    bundle = {"version": "2.0", "snps": snps, "ld": ld}
    with open(out_pkl, "wb") as f:
        pickle.dump(bundle, f, protocol=pickle.HIGHEST_PROTOCOL)

    # quick quality readout
    n_opts = sum(len(s["options"]) for s in snps.values())
    n_gwas = sum(len(s["gwas"]) for s in snps.values())
    n_cat = sum(len(s["catalog"]) for s in snps.values())
    print(f"\nreference.pkl written: {out_pkl}")
    print(f"  SNPs ............. {len(snps):,}")
    print(f"  genotype options . {n_opts:,}")
    print(f"  SNP-page GWAS .... {n_gwas:,}")
    print(f"  catalog assocs ... {n_cat:,}")
    return bundle


if __name__ == "__main__":
    args = sys.argv[1:]
    snp = args[0] if len(args) > 0 else "data/snp_df.csv"
    geno = args[1] if len(args) > 1 else "data/geno_df.csv"
    trait = args[2] if len(args) > 2 else "data/trait_df.csv"
    eq = args[3] if len(args) > 3 else "data/equilibrium_df.csv"
    out = args[4] if len(args) > 4 else "data/reference.pkl"
    # optional 6th arg: deana evidence-pack shards dir; 7th: evidence levels for
    # the broad merge ("high,moderate" default; "high" leaner; "none" to disable)
    deana = args[5] if len(args) > 5 else None
    if len(args) > 6:
        lv = args[6].strip().lower()
        deana_levels = () if lv in ("none", "off", "") else tuple(lv.split(","))
    else:
        deana_levels = ("high", "moderate")
    build(snp, geno, trait, eq, out, deana_dir=deana, deana_levels=deana_levels)
