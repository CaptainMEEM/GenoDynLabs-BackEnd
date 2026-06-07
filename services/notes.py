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


# Strip any allele / genotype / odds-ratio jargon out of free-text summaries so
# the reader sees a plain implication, never "the G allele" or "OR 1.31".
_SCRUB_PATTERNS = [
    re.compile(r"(?i)\b(\d+\s+)?cop(y|ies)\s+of\s+the\s+[ACGT]{1,2}\s+(risk\s+)?allele"),
    re.compile(r"(?i)\bthe\s+[ACGT]{1,2}\s+(risk\s+)?allele\b"),
    re.compile(r"(?i)\b[ACGT]{1,2}\s+risk\s+allele\b"),
    re.compile(r"(?i)\b(odds\s+ratio|OR)\s*[:=]?\s*\d+(\.\d+)?"),
    re.compile(r"(?i)\bp\s*[=<]\s*\d[\d.eE+-]*"),
    re.compile(r"(?i)\bgenotyp\w*\b"),
    # strip any supplementation suggestion (user: no supplement notation anywhere)
    re.compile(r"(?i)\s*[;,\u2014\-]?\s*\b(consider\s+)?supplement\w*\b[^.?!]*[.?!]?"),
    re.compile(r"\((?:\s*[\u2191\u2193]\s*)?\)"),     # empty arrow parens
]


def _scrub(text):
    if not text:
        return ""
    s = text
    for rx in _SCRUB_PATTERNS:
        s = rx.sub("", s)
    s = re.sub(r"\s*[\u2014\-]\s*$", "", s)       # dangling dash
    s = re.sub(r"\s+([:;,.])", r"\1", s)          # space before punctuation
    s = re.sub(r"[:;,]\s*$", "", s)               # dangling colon/comma
    s = re.sub(r"\s{2,}", " ", s).strip(" .;,-\u2014")
    return s


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


# Words that mark a "level/biomarker" trait (phrased as higher/lower levels)
# vs. a disease/risk trait (phrased as an X-times risk).
_LEVEL_RX = re.compile(
    r"(?i)\b(level|levels|concentration|status|biomarker|metabolite|"
    r"vitamin|folate|homocysteine|ferritin|cholesterol|amino acid|"
    r"glycosylation|fatty acid)\b")


def _clean_trait(trait):
    """Tidy a trait string for human reading."""
    t = (trait or "").split(" association near")[0].strip()
    t = re.sub(r"\s+", " ", t).strip(" .;,-")
    return t


# State words that already describe a condition -> never append "levels" to them.
_STATE_RX = re.compile(
    r"(?i)(insufficien|deficien|disease|syndrome|intoleran|disorder|cancer|"
    r"\brisk\b|tolerance|persistence|flush|sensitivity)")

# Stop-words that don't identify the topic of a trait (used for de-duplication).
_STOP = {"level", "levels", "status", "concentration", "plasma", "serum",
         "blood", "total", "risk", "trait", "traits", "disease", "association",
         "circulating", "biomarker", "biomarkers", "values", "value", "higher",
         "lower", "with"}


def _core_terms(text):
    """The topic-identifying words in a trait/summary, for de-dup by subject."""
    return {w for w in re.findall(r"[a-z0-9]+", (text or "").lower())
            if len(w) > 3 and w not in _STOP}


def _noun_form(trait):
    """Trait phrased as a noun for 'higher/lower ___'. Adds 'levels' only when
    the trait is a measurable quantity, never to a state word like 'deficiency'."""
    low = trait.lower()
    if _STATE_RX.search(low) or "level" in low or "status" in low or \
       "concentration" in low:
        return trait
    return trait + " levels"


def _risk_multiplier(orv):
    """Turn an odds ratio into a plain 'X times higher/lower risk' phrase, but
    only when it is a believable case/control OR (roughly 1.1-10x). Returns ''
    when the number isn't a trustworthy risk multiplier."""
    if not orv or orv <= 0:
        return ""
    if 1.1 <= orv <= 10:
        return f"about {orv:.1f}\u00d7 higher risk"
    if 0.1 <= orv <= 0.91:
        return f"about {1.0 / orv:.1f}\u00d7 lower risk"
    return ""               # ~1.0 (no effect) or implausible -> no multiplier


def _implication(trait, risk, orv, user_geno, allow_typical=True):
    """Plain-language implication for one association, with NO allele/copy talk.

    - Level/biomarker traits  -> 'Linked to higher/lower <trait>' (direction from
      the odds ratio); 'Typical <trait>' when the user doesn't carry the variant
      (only if allow_typical).
    - Disease/risk traits     -> 'About 1.5x higher risk of <disease>' when the
      user carries the effect allele and the OR is a believable multiplier;
      otherwise 'Associated with <disease>'. Never states a risk for a variant
      the user does not carry.
    Trait case is preserved (so 'LDL', 'Alzheimer' stay correct).
    """
    trait = _clean_trait(trait)
    if not trait:
        return ""
    dose = _dosage(user_geno, risk)            # used only to decide IF it applies
    is_level = bool(_LEVEL_RX.search(trait)) and not _STATE_RX.search(trait)

    if is_level:
        noun = _noun_form(trait)
        if dose is None:                       # no allele info -> generic
            return f"Associated with {noun}"
        if dose == 0:
            return f"Typical {noun}" if allow_typical else ""
        if orv and orv > 1:
            return f"Linked to higher {noun}"
        if orv and orv < 1:
            return f"Linked to lower {noun}"
        return f"Linked to altered {noun}"

    # disease / risk / state trait
    if dose == 0:
        return ""                              # user doesn't carry it -> don't state
    mult = _risk_multiplier(orv)
    if mult and dose:
        return f"{mult[0].upper()}{mult[1:]} of {trait}"
    if dose or dose is None:
        return f"Associated with {trait}"
    return ""


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

    # ClinVar-style pathogenicity carried on any catalog entry (from deana).
    pathogenic = ""
    for a in snp_record.get("catalog", []):
        sig = str(a.get("sig", "") or "")
        low = sig.lower()
        if "pathogenic" in low and "conflict" not in low and "non-patho" not in low:
            pathogenic = sig
            break

    assocs = _collect_assocs(snp_record)
    effect, eff_or, top_trait = _best_effect_allele(assocs, user_geno)

    # ---- assemble the note from the priority layers --------------------
    parts = []

    # Layer 1: a real per-genotype summary leads (most specific, already plain
    # language e.g. "somewhat lower vitamin B12 levels" / "65% efficiency...").
    clean_summary = _scrub(summary)
    have_summary = not _is_boiler(clean_summary)
    if have_summary:
        parts.append(clean_summary.rstrip("."))

    # Track the topics already covered so we never repeat the same subject
    # (e.g. "lower vitamin D" then "vitamin D insufficiency").
    covered = _core_terms(clean_summary)

    # Layer 2/3: trait associations as PLAIN-LANGUAGE implications (no alleles,
    # no copies, no bare odds ratios). Add up to 2 on NEW topics.
    added = 0
    for a in assocs:
        if added >= 2:
            break
        ct = _clean_trait(a["trait"])
        terms = _core_terms(ct)
        if terms and terms & covered:           # same subject already covered
            continue
        sent = _implication(a["trait"], a.get("risk"), a.get("or"), user_geno,
                            allow_typical=not have_summary)
        if not sent:
            continue
        parts.append(sent)
        covered |= terms
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
            parts.append(_scrub(desc))

    # Final fallback: never blank. State what we do know.
    if not parts:
        if top_trait:
            parts.append(f"Variant studied in relation to {_clean_trait(top_trait)}")
        elif gene:
            parts.append(f"Variant in {gene}")
        else:
            parts.append("Variant of unknown significance")

    # Capitalize the first letter of each part, then join.
    parts = [p[0].upper() + p[1:] if p else p for p in parts]
    note = ". ".join(p for p in parts if p).strip()
    if note and not note.endswith((".", "!", "?")):
        note += "."

    # Surface clinical pathogenicity prominently.
    if pathogenic:
        tag = "Likely pathogenic" if "likely" in pathogenic.lower() and \
              "/" not in pathogenic else pathogenic
        note = f"ClinVar: {tag}. " + note

    hl = _highlight(effect, user_geno, repute)
    if pathogenic and user_geno not in ("--", ""):
        hl = "hl-orange"

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
        "pathogenic": bool(pathogenic),
        "label": gene_labels.label(gene, fallback_trait=top_trait),
        "label_known": gene_labels.is_known(gene),
    }
