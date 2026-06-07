"""
notes.py  --  Genotype-FIRST note synthesis.

Two hard rules drive this module:

  1. CARRIED-ONLY.  A variant is only "informative" for a user when the user's
     own genotype shows something: they carry the effect/risk allele, OR SNPedia
     has a real (non-boiler) summary for *their* exact genotype, OR the genotype
     is clinically/SNPedia-flagged. Reference/wild-type rows with no signal are
     marked `informative=False` and the report drops them. This is the fix for
     "we show rows the user doesn't even have" — the Effect allele we display is
     ALWAYS one the user actually carries (else we show no effect allele at all),
     so "Effect: G / Your Genotype: CC" can no longer happen.

  2. ACTIONABLE NOTES.  The note describes what the user's genotype *means*, not
     what the gene does in the abstract (that context is the article subtitle).
     Carrying a nutrient allele reads like "One copy of the A allele — associated
     with slightly lower Vitamin B6 levels", not "Part of the pathway that ...".

Public entry point:
    annotate_variant(user_geno, snp_record, genotype_option) -> dict

Strand note: `user_geno` arrives already strand-corrected to SNPedia orientation
and sorted (e.g. "AG"). SNPedia genotype summaries are in that same orientation,
so the highest-quality path is reliable. For a GWAS *risk allele* we only ever
claim carriage when that allele literally appears in the normalized genotype;
strand-ambiguous cases are treated as "not carried" rather than guessed, which
is what keeps nonsense effect/genotype pairs off the report.
"""
import re

try:
    from . import gene_labels
except ImportError:                       # standalone / testing
    import gene_labels

ARROW_UP = "\u2191"     # ↑
ARROW_DN = "\u2193"     # ↓
EMDASH = "\u2014"
_VALID = set("ACGT")

# Genotype summaries that carry no biological signal. Matched as the WHOLE
# string so a real note that merely begins with "Common variant..." survives.
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


def _dosage(user_geno, allele):
    """Copies of `allele` the user carries (0/1/2), or None if unknown.
    Only counts when `allele` is a single valid base actually present-countable."""
    if (not allele) or (not user_geno) or user_geno in ("--", ""):
        return None
    if len(allele) != 1 or allele not in _VALID:
        return None
    return user_geno.count(allele)


def _direction_word(orv):
    """'higher' / 'lower' / 'slightly higher' / 'slightly lower' / '' (unknown)."""
    if not orv or orv <= 0:
        return ""
    if 0.83 <= orv <= 1.20:                 # close to null -> small effect
        return "slightly higher" if orv >= 1 else "slightly lower"
    return "higher" if orv > 1 else "lower"


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
                merged[key] = {"trait": trait, "risk": (a.get("risk") or "").upper(),
                               "or": orv}
            else:
                if keep.get("or") is None and orv is not None:
                    keep["or"] = orv
                if not keep.get("risk") and a.get("risk"):
                    keep["risk"] = a.get("risk").upper()
    return sorted(merged.values(),
                  key=lambda a: abs((a["or"] or 1) - 1), reverse=True)


def _carried_assoc(assocs, user_geno):
    """Return the strongest association whose risk allele the user ACTUALLY
    carries, as (risk, or, trait, dose). If none is carried, return the strongest
    association's risk allele with dose 0 (used only for context, never shown)."""
    best_uncarried = None
    for a in assocs:
        risk = a.get("risk") or ""
        dose = _dosage(user_geno, risk)
        if risk and dose:                   # dose >= 1 -> user carries it
            return risk, a.get("or"), a["trait"], dose
        if best_uncarried is None:
            best_uncarried = (risk, a.get("or"), a.get("trait", ""), dose or 0)
    if best_uncarried:
        return best_uncarried
    return "", None, "", 0


def _dose_phrase(dose, risk):
    if dose >= 2:
        return f"Two copies of the {risk} allele"
    return f"One copy of the {risk} allele"


def _anchor_phrase(nutrient, trait):
    """The thing the genotype affects, phrased for a reader."""
    if nutrient:
        n = nutrient.strip()
        # nutrient labels are like "Vitamin B6", "Iron", "Vitamin C"
        return f"{n} levels" if not n.lower().endswith(("levels", "status")) else n
    t = (trait or "").strip()
    if not t or t.lower() in ("nan", "none", "nr"):
        return ""
    return t[:1].lower() + t[1:]            # keep gene-style traits lowercase-led


def _carrier_sentence(dose, risk, orv, nutrient, trait):
    """A user-specific, directional sentence, e.g.
    'One copy of the A allele — associated with slightly lower Vitamin B6 levels.'"""
    anchor = _anchor_phrase(nutrient, trait)
    head = _dose_phrase(dose, risk)
    direction = _direction_word(orv)
    if anchor and direction:
        return f"{head} {EMDASH} associated with {direction} {anchor}"
    if anchor:
        return f"{head} {EMDASH} linked to {anchor} in genetic studies"
    return f"{head} {EMDASH} flagged in genetic studies"


def _highlight(carried, dose, repute, pathogenic_carried):
    if pathogenic_carried:
        return "hl-orange"
    if carried:
        if dose >= 2:
            return "hl-orange"
        return "hl-yellow"
    rl = (repute or "").strip().lower()
    if rl == "good":
        return "hl-good"
    if rl == "bad":
        return "hl-bad"
    return ""


def annotate_variant(user_geno, snp_record, option):
    """Build the full annotation for one matched variant.

    Returns a dict the report renders directly. Key additions vs. the old
    pipeline: `informative` (drop-if-False), `carried`, and an `effect_allele`
    that is non-empty ONLY when the user carries it.
    """
    gene = snp_record.get("gene", "")
    nutrient = snp_record.get("nutrient")
    summary = (option.get("summary") or "").strip()
    repute = (option.get("repute") or "neutral").strip() or "neutral"
    mag = option.get("mag") or 0.0
    geno_known = bool(user_geno) and user_geno not in ("--", "")

    # ClinVar-style pathogenicity carried on any catalog entry (from deana).
    pathogenic = ""
    for a in snp_record.get("catalog", []):
        sig = str(a.get("sig", "") or "")
        low = sig.lower()
        if "pathogenic" in low and "conflict" not in low and "non-patho" not in low:
            pathogenic = sig
            break

    assocs = _collect_assocs(snp_record)
    risk, eff_or, top_trait, dose = _carried_assoc(assocs, user_geno)
    carried = bool(risk) and dose >= 1

    real_summary = bool(summary) and not _is_boiler(summary)
    flagged = (mag >= 2.0) or ((repute.lower() in ("bad", "good")) and mag >= 1.0)
    pathogenic_carried = bool(pathogenic) and carried

    # The variant earns a row ONLY if the user's own genotype shows something.
    informative = bool(geno_known and (
        carried or real_summary or pathogenic_carried or flagged))

    # Effect allele to DISPLAY: only ever an allele the user actually carries.
    effect_display = risk if carried else ""

    # ---- assemble the note: always about the user's genotype ----------------
    parts = []
    if pathogenic_carried:
        tag = ("Likely pathogenic" if "likely" in pathogenic.lower()
               and "/" not in pathogenic else pathogenic)
        anchor = _anchor_phrase(nutrient, top_trait)
        carrier = (f"you carry {('two copies' if dose >= 2 else 'one copy')}"
                   if carried else "present")
        tail = f" ({anchor})" if anchor else ""
        parts.append(f"ClinVar: {tag}{tail} {EMDASH} {carrier}")

    if real_summary:
        # SNPedia's per-genotype sentence is already specific to THIS genotype.
        parts.append(summary.rstrip("."))
    elif carried:
        parts.append(_carrier_sentence(dose, risk, eff_or, nutrient, top_trait))
    elif flagged and not parts:
        desc = (option.get("desc") or "").strip()
        if desc and not _is_boiler(desc):
            parts.append(desc.rstrip("."))
        else:
            parts.append(f"SNPedia flags your genotype as notable "
                         f"(magnitude {mag:.1f})")

    # Safety net so the field is never blank for a row we *do* show. (Filtered
    # rows keep whatever is here but are not rendered.)
    if not parts:
        if real_summary:
            parts.append(summary.rstrip("."))
        elif geno_known:
            parts.append(f"You carry the {user_geno} genotype")
        else:
            parts.append("Not genotyped in your data")

    note = ". ".join(p for p in parts if p).strip()
    if note and not note.endswith((".", "!", "?")):
        note += "."

    hl = _highlight(carried, dose, repute, pathogenic_carried)

    return {
        "note": note,
        "effect_allele": effect_display,     # carried-only; "" otherwise
        "risk_allele": risk,                 # biological risk allele (for scoring)
        "carried": carried,
        "dose": dose if dose is not None else 0,
        "primary_trait": top_trait,
        "max_or": eff_or if carried else None,
        "mag": mag,
        "repute": repute,
        "hl": hl,
        "gene": gene,
        "pathogenic": bool(pathogenic_carried),
        "informative": informative,
        "label": gene_labels.label(gene, fallback_trait=top_trait),
        "label_known": gene_labels.is_known(gene),
    }
