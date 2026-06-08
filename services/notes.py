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


def _strongest_assoc(assocs, allele):
    """The (or, trait) of the strongest association naming `allele` as its risk
    allele, else the strongest association's (or, trait), else (None, "")."""
    for a in assocs:
        if (a.get("risk") or "").upper() == allele:
            return a.get("or"), a.get("trait", "")
    if assocs:
        return assocs[0].get("or"), assocs[0].get("trait", "")
    return None, ""


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
    rl = (repute or "").strip().lower()
    if pathogenic_carried:
        return "hl-orange"
    if carried:
        if rl == "bad":
            return "hl-bad"                  # carried + adverse keeps its tint
        if rl == "good":
            return "hl-good"
        return "hl-orange" if dose >= 2 else "hl-yellow"
    if rl == "bad":
        return "hl-bad"
    if rl == "good":
        return "hl-good"
    return ""


def _scored_options(snp_record):
    """[(sorted_geno, notability)] for each documented biallelic genotype, where
    notability ranks how 'interesting' SNPedia considers that genotype (magnitude,
    with a non-neutral repute as a secondary signal so flat-magnitude SNPs still
    rank)."""
    out = []
    for geno, opt in (snp_record.get("options") or {}).items():
        g = "".join(ch for ch in (geno or "").upper() if ch in _VALID)
        if len(g) != 2:
            continue
        mag = float(opt.get("mag") or 0.0)
        rep = (opt.get("repute") or "").strip().lower()
        rep_bonus = 0.5 if rep in ("bad", "good") else 0.0
        out.append((g, mag + rep_bonus))
    return out


def _effect_from_options(snp_record):
    """Derive the effect/variant allele from SNPedia's own genotype table: the
    allele that distinguishes the most-notable documented genotype from the
    least-notable (reference) one. Returns "" when it can't be told apart."""
    scored = _scored_options(snp_record)
    if len(scored) < 2:
        return ""
    scored.sort(key=lambda x: x[1])
    lo_geno, lo_n = scored[0]
    hi_geno, hi_n = scored[-1]
    if hi_n <= lo_n:                          # no spread -> nothing to learn
        return ""
    distinguishing = set(hi_geno) - set(lo_geno)
    return next(iter(distinguishing)) if len(distinguishing) == 1 else ""


def _is_reference_geno(snp_record, user_geno):
    """True only when we can confidently say the user's genotype is the SNP's
    reference (least-notable) genotype. False if the options can't tell us."""
    scored = _scored_options(snp_record)
    if len(scored) < 2:
        return False
    scored.sort(key=lambda x: x[1])
    return scored[0][0] == "".join(sorted(c for c in user_geno if c in _VALID))


def _variant_allele(snp_record, assocs):
    """The SNP's single effect/variant allele, from the strongest available
    signal: (1) SNPedia magnitude/repute table, then (2) any GWAS/Catalog risk
    allele. Independent of the user — carriage is checked separately. "" if the
    SNP gives us no way to name a variant allele."""
    ea = _effect_from_options(snp_record)
    if ea:
        return ea
    for a in assocs:
        r = (a.get("risk") or "").upper()
        if len(r) == 1 and r in _VALID:
            return r
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
    real_summary = bool(summary) and not _is_boiler(summary)

    # UNIVERSAL EFFECT-ALLELE RULE.
    # Identify the SNP's variant (effect) allele from any signal, then keep the
    # row only if the user CARRIES it. The Effect column is therefore always a
    # single allele that is present in the user's genotype — never blank on a
    # shown row, and never an allele the user doesn't have.
    risk = _variant_allele(snp_record, assocs)

    # Last-resort: a documented finding on a HOMOZYGOUS genotype means both
    # copies are the relevant allele, so that allele is the effect — unless the
    # options table tells us this genotype is actually the reference one.
    if not risk and geno_known:
        bases = [c for c in user_geno if c in _VALID]
        if len(bases) == 2 and bases[0] == bases[1] and (
                real_summary or pathogenic or repute.lower() in ("bad", "good")):
            if not _is_reference_geno(snp_record, user_geno):
                risk = bases[0]

    dose = _dosage(user_geno, risk) or 0
    carried = bool(risk) and dose >= 1
    eff_or, top_trait = _strongest_assoc(assocs, risk) if carried else (None, "")
    pathogenic_carried = bool(pathogenic) and carried

    # A row is shown ONLY when we can name an effect allele the user carries.
    # No carried effect allele -> no row (this is "if it doesn't match the
    # genotype it shouldn't be there", applied uniformly to every variant).
    informative = bool(geno_known and carried)

    # Effect allele to DISPLAY: always the carried allele on a shown row.
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
        # SNPedia's per-genotype sentence is already specific to THIS genotype;
        # the Effect column (now filled above) shows which allele it concerns.
        parts.append(summary.rstrip("."))
    elif carried:
        parts.append(_carrier_sentence(dose, risk, eff_or, nutrient, top_trait))

    # Safety net so the field is never blank for a row we *do* show. (Filtered
    # rows keep whatever is here but are not rendered.)
    if not parts:
        if geno_known:
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
