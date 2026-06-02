"""
Run ONCE at build/deploy time (not per request).

Parses the giant SNPedia dumps (snp_df.csv ~33MB of regex-laden text, plus
geno_df.csv) into a compact rsid-keyed dict and pickles it. After this, no
request ever touches the 33MB raw files or runs the parsing regexes again.

    python precompute_reference.py data/snp_df.csv data/geno_df.csv data/reference.pkl
"""
import re, sys, pickle
import pandas as pd


def field(name, text):
    m = re.search(r"%s=([^$|}]+)" % re.escape(name), text or "")
    return m.group(1).strip() if m else ""


def clean(text, limit=180):
    s = re.sub(r"\{\{[^}]*\}\}", " ", text or "")
    s = re.sub(r"\[\[([^\]]+)\]\]", r"\1", s)
    s = s.replace("$", " ")
    s = re.sub(r"\s+", " ", s).strip(" .;,-")
    return s[:limit]


def build(snp_csv, geno_csv, out_pkl):
    # SNP-level info (gene/chr/pos/orientation), keyed by rsid
    snp = {}
    for t in pd.read_csv(snp_csv)["text"].dropna():
        rs = field("rsid", t)
        if not rs:
            continue
        snp["rs" + rs] = {
            "gene": field("Gene", t),
            "chr": field("Chromosome", t),
            "pos": field("position", t),
            "orientation": (field("StabilizedOrientation", t)
                            or field("Orientation", t)).lower(),
        }

    def parse_opt(text):
        a1, a2 = field("allele1", text), field("allele2", text)
        geno = "".join(sorted([a1, a2])).lstrip("=")
        summ = field("summary", text) or clean(text)
        return {"geno": geno, "mag": field("magnitude", text) or "0",
                "repute": field("repute", text) or "neutral",
                "summary": summ.replace("\n", " "), "rsid": field("rsid", text)}

    ref = {}
    for rec in pd.read_csv(geno_csv).itertuples():
        opts = [parse_opt(str(getattr(rec, c))) for c in ("gt1", "gt2", "gt3")]
        rs = next(("rs" + o["rsid"] for o in opts if o["rsid"]), "")
        if not rs:
            continue
        base = snp.get(rs, {"gene": "", "chr": "", "pos": "", "orientation": ""})
        ref[rs] = {
            "gene": base["gene"], "chr": base["chr"], "pos": base["pos"],
            "orientation": base["orientation"],
            "options": [{"geno": o["geno"], "mag": o["mag"],
                         "repute": o["repute"], "summary": o["summary"]}
                        for o in opts if o["geno"]],
        }

    with open(out_pkl, "wb") as f:
        pickle.dump(ref, f, protocol=pickle.HIGHEST_PROTOCOL)
    return ref


if __name__ == "__main__":
    snp = sys.argv[1] if len(sys.argv) > 1 else "/mnt/user-data/uploads/snp_df.csv"
    geno = sys.argv[2] if len(sys.argv) > 2 else "/mnt/user-data/uploads/geno_df.csv"
    out = sys.argv[3] if len(sys.argv) > 3 else "reference.pkl"
    r = build(snp, geno, out)
    print(f"reference: {len(r):,} rsids -> {out}")
