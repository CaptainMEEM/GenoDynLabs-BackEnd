"""
In-process genome annotator. Loads the precomputed reference ONCE per worker
(module-level cache) and matches a user's 23andMe file with O(1) dict lookups.
Replaces the per-request subprocess + O(N^2) scan in generate_variant_report.py.

    from annotate import annotate_genome
    n = annotate_genome(clean_path, out_csv="snpedia_data.csv",
                        ref_path="data/reference.pkl")
"""
import csv, pickle

_REF = None
_COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}
_COLS = ["rsid", "geno", "mag", "repute", "summary",
         "description", "geno_description", "chr", "pos", "gene"]


def load_reference(ref_path="data/reference.pkl"):
    """Load once; reused for every subsequent request in this process."""
    global _REF
    if _REF is None:
        with open(ref_path, "rb") as f:
            _REF = pickle.load(f)
    return _REF


def _flip(g):
    return "".join(_COMP.get(b, b) for b in g)[::-1]


def annotate_genome(genome_path, out_csv, ref_path="data/reference.pkl"):
    ref = load_reference(ref_path)

    # Build user rsid->genotype only for rsids SNPedia knows (intersection).
    user = {}
    with open(genome_path, newline="", errors="replace") as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r, None) or []
        try:
            ri, gi = header.index("rsid"), header.index("genotype")
        except ValueError:
            ri, gi = 0, 3
        wide = max(ri, gi)
        for row in r:
            if len(row) <= wide:
                continue
            rsid = row[ri].strip()
            if rsid in ref:                      # O(1)
                user[rsid] = row[gi].strip().upper()

    n = 0
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=_COLS)
        w.writeheader()
        for rsid, ug in user.items():
            rec = ref[rsid]                      # O(1)
            g = _flip(ug) if rec["orientation"] == "minus" else ug
            g = "".join(sorted(g))
            for o in rec["options"]:
                if o["geno"] == g:
                    w.writerow({"rsid": rsid, "geno": g, "mag": o["mag"],
                                "repute": o["repute"], "summary": o["summary"],
                                "description": "", "geno_description": "",
                                "chr": rec["chr"], "pos": rec["pos"],
                                "gene": rec["gene"]})
                    n += 1
                    break
    return n
