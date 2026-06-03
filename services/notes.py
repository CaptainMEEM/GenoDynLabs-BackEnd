"""
notes.py  --  Layered note synthesis. The engine that fixes "bad notes".

For every matched variant we assemble the single best note we can from ALL of
SNPedia/GWAS, in priority order, and we GUARANTEE a meaningful note for every
gene pair. The layers, highest priority first:

  1. Genotype summary      (geno_df)  -- the per-genotype sentence, when it is
                                          real (not "common in clinvar" filler).
  2. SNP-page GWAS          (snp_df)  -- named trait + your-dosage direction + OR.
  3. GWAS-Catalog assoc.  (trait_df)  -- study-backed trait + risk allele + OR.
  4. Curated gene function (gene_labels) -- always-true "what this gene does".
  5. Genotype-page desc.   (geno_df)  -- cleaned page text, last resort.

We also compute the effect allele to show, the user's dosage of it, and a
highlight class, all from the same evidence the note is built from -- so the
"Effect" / "Your Genotype" columns and the note always agree.

Public entry point:
    annotate_variant(user_geno, snp_record, genotype_option) -> dict
"""
import re

try:
    from . import gene_labels
except ImportError:                       # standalone / testing
    import gene_labels

ARROW_UP = "\u2191"     # ↑
ARROW_DN = "\u2193"     # ↓

# Genotype summaries that carry no biological signal. Matched as the WHOLE
# string so a real note that merely begins with "Common variant..." survives.
# "common in <anything>" (clinvar / complete genomics / 1000 genomes / ...) and
# bare "normal"/"benign"/"none" are all filler.
_BOILER = re.compile(
    r"(?i)^\s*("
    r"common(\s+(in|on)\b.*|/\w+)?|"
    r"normal|none|n/?a|unknown|"
    r"benign|likely\s+benign|uncertain\s+significance|"
    r"no\s+\w+|not\s+\w+|wild\s*type|reference"
    r")\s*\.?\s*$"
)


def _is_boiler(summary):
    return (not summary) or bool(_BOILER.match(summary))


def _or_phrase(orv):
    """Readable OR clause. Values outside [0.1, 10] are almost always effect
    sizes on continuous traits (metabolite/biomarker levels), not true
    case/control odds ratios, so we show direction only rather than a
    misleading number like 'OR 49.77'. (Same bounds ModernPromethease used.)"""
    if not orv:
        return ""
    arrow = ARROW_DN if orv < 1 else ARROW_UP
    if orv < 0.1 or orv > 10:
        return arrow            # direction only; magnitude not trustworthy
    return f"{arrow} OR {orv:.2f}"


def _dosage(user_geno, allele):
    """How many copies of `allele` the user carries (0/1/2), or None if unknown."""
    if not allele or not user_geno or user_geno in ("--", ""):
        return None
    return user_geno.count(allele)


def _assoc_sentence(trait, risk, orv, user_geno):
    """One readable clause for a single trait association, dosage-aware."""
    trait = (trait or "").strip()
    if not trait:
        return ""
    dose = _dosage(user_geno, risk)
    pieces = []
    if dose is not None and risk:
        if dose == 0:
            pieces.append(f"you carry no {risk} risk allele")
        elif dose == 1:
            pieces.append(f"1 copy of the {risk} risk allele")
        else:
            pieces.append(f"2 copies of the {risk} risk allele")
    orp = _or_phrase(orv)
    head = trait
    if orp:
        head = f"{trait} ({orp})"
    if pieces:
        return f"{head} \u2014 {pieces[0]}"
    return head


def _collect_assocs(snp_record):
    """Merge SNP-page GWAS + GWAS-Catalog associations, strongest first,
    de-duplicated by trait. Catalog (study-backed) wins ties."""
    merged = {}
    for src in (snp_record.get("catalog", []), snp_record.get("gwas", [])):
        for a in src:
            trait = (a.get("trait") or "").strip()
            if not trait:
                continue
            key = trait.lower()
            orv = a.get("or")
            keep = merged.get(key)
            if keep is None:
                merged[key] = {"trait": trait, "risk": a.get("risk", ""), "or": orv}
            else:
                # prefer an entry that actually has an OR / risk allele
                if keep.get("or") is None and orv is not None:
                    keep["or"] = orv
                if not keep.get("risk") and a.get("risk"):
                    keep["risk"] = a.get("risk")
    return sorted(merged.values(),
                  key=lambda a: abs((a["or"] or 1) - 1), reverse=True)


def _best_effect_allele(assocs, user_geno):
    """Pick the effect allele to display: the strongest association's risk
    allele that the SNP actually carries information about."""
    for a in assocs:
        if a.get("risk"):
            return a["risk"], a.get("or"), a["trait"]
    return "", None, (assocs[0]["trait"] if assocs else "")


def _highlight(effect_allele, user_geno, repute):
    if user_geno in ("--", ""):
        return ""
    if effect_allele:
        dose = user_geno.count(effect_allele)
        if dose >= 2:
            return "hl-orange"
        if dose == 1:
            return "hl-yellow"
        return ""
    if (repute or "").strip().lower() == "bad":
        return "hl-orange"
    return ""


def annotate_variant(user_geno, snp_record, option):
    """Build the full annotation for one matched variant.

    user_geno   : the user's genotype, strand-corrected & sorted, e.g. "AG"
    snp_record  : the rich per-SNP dict from reference.pkl
    option      : the matched genotype option dict (mag/repute/summary/desc)

    Returns a dict the report renders directly.
    """
    gene = snp_record.get("gene", "")
    summary = (option.get("summary") or "").strip()
    repute = (option.get("repute") or "neutral").strip() or "neutral"
    mag = option.get("mag") or 0.0

    assocs = _collect_assocs(snp_record)
    effect, eff_or, top_trait = _best_effect_allele(assocs, user_geno)

    # ---- assemble the note from the priority layers --------------------
    parts = []

    # Layer 1: a real per-genotype summary leads (most specific to the user).
    if not _is_boiler(summary):
        parts.append(summary.rstrip("."))

    # Layer 2/3: trait associations with dosage + direction. Add up to 2 that
    # aren't already implied by the summary.
    added = 0
    summary_low = summary.lower()
    for a in assocs:
        if added >= 2:
            break
        if a["trait"].lower() in summary_low:
            continue
        sent = _assoc_sentence(a["trait"], a.get("risk"), a.get("or"), user_geno)
        if sent:
            parts.append(sent)
            added += 1

    # Layer 4: curated baseline gene function -- guarantees a note for known
    # genes even when nothing genotype-specific exists.
    if not parts:
        bn = gene_labels.baseline_note(gene)
        if bn:
            parts.append(bn)

    # Layer 5: cleaned genotype-page description, true last resort.
    if not parts:
        desc = (option.get("desc") or "").strip()
        if desc and not _is_boiler(desc):
            parts.append(desc)

    # Final fallback: never blank. State what we do know.
    if not parts:
        if top_trait:
            parts.append(f"Variant studied in relation to {top_trait.lower()}")
        elif gene:
            parts.append(f"Variant in {gene}")
        else:
            parts.append("Variant of unknown significance")

    note = ". ".join(p for p in parts if p).strip()
    if note and not note.endswith((".", "!", "?")):
        note += "."

    hl = _highlight(effect, user_geno, repute)

    return {
        "note": note,
        "effect_allele": effect,
        "dose": _dosage(user_geno, effect) or 0,
        "primary_trait": top_trait,
        "max_or": eff_or,
        "mag": mag,
        "repute": repute,
        "hl": hl,
        "gene": gene,
        "label": gene_labels.label(gene, fallback_trait=top_trait),
        "label_known": gene_labels.is_known(gene),
    }
