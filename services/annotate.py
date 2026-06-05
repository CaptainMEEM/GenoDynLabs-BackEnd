"""
annotate.py  --  In-process genome annotator (rich edition).

Loads the rich reference bundle ONCE per worker (module cache) and matches a
user's 23andMe file with O(1) dict lookups, producing fully-annotated records
(note, effect allele, dosage, highlight, topic-ready label) via notes.py.

    from services.annotate import annotate_genome
    records, stats = annotate_genome(genome_text, ref_path="data/reference.pkl")

`records` is a list of dicts ready for the report; nothing is written to disk.
"""
import pickle

try:
    from . import notes
except ImportError:
    import notes

_BUNDLE = None
_COMP = {"A": "T", "T": "A", "G": "C", "C": "G", "-": "-"}


def load_reference(ref_path="data/reference.pkl"):
    """Load once; reused for every subsequent request in this process."""
    global _BUNDLE
    if _BUNDLE is None:
        with open(ref_path, "rb") as f:
            _BUNDLE = pickle.load(f)
    return _BUNDLE


def _flip(g):
    return "".join(_COMP.get(b, b) for b in g)[::-1]


def _norm(g):
    g = (g or "").strip().upper()
    return "".join(sorted(g)) if g else ""


def _iter_genome(genome_text):
    """Yield (rsid, genotype) from a cleaned 23andMe text blob (tab-separated)."""
    lines = genome_text.splitlines()
    if not lines:
        return
    header = lines[0].lower().split("\t")
    try:
        ri = header.index("rsid")
        gi = header.index("genotype")
    except ValueError:
        ri, gi = 0, 3
    wide = max(ri, gi)
    for line in lines[1:]:
        if not line:
            continue
        row = line.split("\t")
        if len(row) <= wide:
            continue
        yield row[ri].strip(), row[gi].strip().upper()


def annotate_genome(genome_text, ref_path="data/reference.pkl"):
    """Match the user's genome against the reference and annotate every hit.

    Returns (records, stats). `records` is a list of dicts; each carries the
    raw match plus everything notes.annotate_variant produced.
    """
    bundle = load_reference(ref_path)
    snps = bundle["snps"]

    records = []
    n_seen = n_matched = 0
    for rsid, raw_geno in _iter_genome(genome_text):
        rec = snps.get(rsid) or snps.get(rsid.lower())
        if rec is None:
            continue
        n_seen += 1

        # strand-correct against SNPedia orientation, then sort alleles
        ug = _flip(raw_geno) if rec.get("orientation") == "minus" else raw_geno
        ug = _norm(ug)

        # Find the matching genotype option. If SNPedia has no exact option for
        # this genotype we still keep the variant (with a neutral option) so the
        # SNP's GWAS/curated note can carry it -- this is a big source of the
        # coverage the old pipeline lost.
        option = rec["options"].get(ug)
        if option is None:
            # try the raw (un-flipped) orientation as a fallback match
            option = rec["options"].get(_norm(raw_geno))
        if option is None:
            option = {"mag": 0.0, "repute": "neutral", "summary": "", "desc": ""}

        ann = notes.annotate_variant(ug or _norm(raw_geno), rec, option)
        ann.update({
            "rsid": rsid,
            "geno": ug or _norm(raw_geno) or "--",
            "chr": rec.get("chr", ""),
            "pos": rec.get("pos", ""),
            "genes": rec.get("genes", []),
            "nutrient": rec.get("nutrient"),     # deana-tagged panel routing
        })
        records.append(ann)
        n_matched += 1

    stats = {"genotyped_in_snpedia": n_seen, "annotated": n_matched}
    return records, stats
