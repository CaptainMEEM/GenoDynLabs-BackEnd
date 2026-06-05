"""
features.py  --  Curated "features": named diet / methylation / intolerance /
condition articles, on top of the nutrient panel.

Each feature groups a set of genes under one article in a chosen topic, e.g.
"Folate Conversion (MTHFR)", "Lactose Tolerance", "CoQ10 Deficiency". A feature
article appears only when the user actually matches one of its genes (so the
report stays personal), except the nutrient panel which is pinned elsewhere.

A gene belongs to at most ONE feature (first match wins), and feature routing
takes PRECEDENCE over the nutrient-panel grouping. That is what restores the
prominent "MTHFR / folate / methylation" framing: MTHFR/MTR/MTRR/CBS leave the
vitamin panel and get their own Methylation section, while the B9/B12 panel
entries keep the absorption/transport genes.

Two archetypes (Mediterranean diet, Hunter-gatherer vs Farmer) are handled
specially in report.py: they are matched-only SYNTHESIZED articles keyed to a
few well-known marker SNPs, because their genes are "owned" by other features.

Each feature: (article_label, topic_key, [genes], baseline_note).
Add freely; genes are upper-cased.
"""

FEATURES = [
    # ---- Methylation & Homocysteine (the famous MTHFR/folate story) --------
    ("Folate Conversion (MTHFR)", "methylation", ["MTHFR"],
     "MTHFR converts folate (folic acid) to its active 5-MTHF form; the C677T/"
     "A1298C variants slow this and can raise homocysteine."),
    ("B12 Methylation (MTR & MTRR)", "methylation", ["MTR", "MTRR"],
     "MTR and MTRR use vitamin B12 to recycle homocysteine back to methionine, "
     "driving the methylation cycle."),
    ("Homocysteine Clearance (CBS)", "methylation", ["CBS", "BHMT"],
     "CBS clears homocysteine down the transsulfuration pathway; BHMT offers a "
     "betaine-dependent backup route."),

    # ---- Diet, Food & Intolerances ----------------------------------------
    ("Lactose Tolerance", "diet", ["MCM6", "LCT"],
     "The MCM6/LCT region controls whether the lactase enzyme stays active in "
     "adulthood, i.e. whether you digest milk sugar."),
    ("Caffeine Metabolism", "diet", ["CYP1A2", "AHR"],
     "CYP1A2 (regulated by AHR) clears most caffeine; variants make people fast "
     "or slow caffeine metabolizers."),
    ("Alcohol Metabolism & Flush", "diet", ["ALDH2", "ADH1B", "ADH1C"],
     "ADH1B/ADH1C set how fast alcohol becomes acetaldehyde and ALDH2 clears it; "
     "the ALDH2 variant causes the alcohol flush. Note: the same acetaldehyde "
     "and the histamine in many drinks can also worsen histamine intolerance."),
    ("Histamine Intolerance", "diet", ["AOC1", "HNMT", "ABP1"],
     "AOC1 (DAO) breaks down histamine from food in the gut and HNMT clears it "
     "inside cells; low activity can cause histamine-intolerance symptoms."),
    ("Tyramine Sensitivity", "diet", ["MAOA", "MAOB"],
     "MAOA/MAOB break down tyramine from aged and fermented foods; lower "
     "activity can mean more sensitivity to tyramine (e.g. headaches)."),
    ("Oxalate Handling", "diet", ["AGXT", "GRHPR", "HOGA1"],
     "These enzymes keep the body's oxalate low; loss-of-function variants raise "
     "oxalate and the risk of calcium-oxalate kidney stones."),
    ("Saturated Fat Response", "diet", ["APOA2"],
     "The APOA2 -265 variant links a high saturated-fat intake to greater weight "
     "gain in carriers."),
    ("Omega-3 / Fish Oil Response", "diet", ["FADS1", "FADS2", "ELOVL2"],
     "FADS1/FADS2 build long-chain omega-3 (EPA/DHA) from plant precursors; "
     "low-efficiency variants benefit more from preformed fish-oil omega-3."),
    ("Amylase / Starch Digestion", "diet", ["AMY1A", "AMY1B", "AMY1C", "AMY2A", "AMY2B"],
     "The amylase region sets how well salivary/pancreatic amylase breaks down "
     "dietary starch (note: the main driver is AMY1 copy number, not single "
     "SNPs, so single-variant results here are only a partial readout)."),
    ("Mushroom / Hydrazine Detox", "diet", ["NAT2"],
     "Acetylator status (NAT2) affects how the body processes various dietary "
     "and environmental amines, including naturally occurring mushroom "
     "compounds. Interpretive only \u2014 not an established 'mushroom' result."),

    # ---- Inherited / metabolic conditions ---------------------------------
    ("Fructose Intolerance", "conditions", ["ALDOB", "KHK", "SLC2A5"],
     "ALDOB breaks down fructose; pathogenic variants cause hereditary fructose "
     "intolerance, where fructose/sucrose must be avoided."),
    ("Creatine Synthesis", "conditions", ["GATM", "GAMT", "SLC6A8"],
     "GATM/GAMT make creatine and SLC6A8 transports it; defects cause cerebral "
     "creatine-deficiency syndromes."),
    ("CoQ10 (Ubiquinone) Deficiency", "conditions",
     ["COQ2", "COQ4", "COQ6", "COQ7", "COQ8A", "COQ8B", "COQ9",
      "PDSS1", "PDSS2", "NQO1", "ADCK3"],
     "These genes synthesize coenzyme Q10; variants can lower CoQ10 and may "
     "inform whether CoQ10 supplementation is worth discussing."),
    ("Wilson's Disease (Copper)", "conditions", ["ATP7B"],
     "ATP7B exports excess copper; pathogenic variants cause Wilson's disease, "
     "a treatable copper-overload disorder."),
    ("Hunter Syndrome (MPS II)", "conditions", ["IDS"],
     "IDS breaks down certain complex sugars; pathogenic variants cause Hunter "
     "syndrome (mucopolysaccharidosis type II)."),

    # ---- Allergy / atopy (immune section) ---------------------------------
    ("Allergy & Atopy", "immune", ["FLG", "IL13", "IL4", "IL4R", "TSLP", "IL33"],
     "These barrier/immune genes influence eczema, allergic sensitization, and "
     "atopic conditions."),
    ("Nickel Allergy", "immune", ["TNF"],
     "Immune variants (e.g. TNF, HLA) associated with susceptibility to contact "
     "allergy, including nickel contact dermatitis. Interpretive."),
]

# Archetype panels: matched-only, synthesized from specific marker SNPs in
# report.py. marker rsid -> short descriptor.
ARCHETYPES = {
    "Mediterranean Diet Response": {
        "topic": "diet",
        "note": "These variants modify how much someone benefits from a "
                "Mediterranean-style diet (olive oil, fish, vegetables).",
        "markers": {
            "rs7903146": "TCF7L2 \u2014 carriers gain more cardiometabolic "
                         "benefit from a Mediterranean diet",
            "rs5082": "APOA2 \u2014 modifies the saturated-fat/Mediterranean "
                      "response",
            "rs9939609": "FTO \u2014 Mediterranean diet blunts this weight-gain "
                         "variant",
            "rs1799983": "NOS3 \u2014 interacts with a Mediterranean diet for "
                         "blood-pressure benefit",
            "rs5186": "AGTR1 \u2014 modifies blood-pressure response to a "
                      "Mediterranean diet",
        },
    },
    "Hunter-Gatherer vs Farmer (interpretive)": {
        "topic": "diet",
        "note": "An interpretive archetype: variants in starch, dairy, and fat "
                "adaptation that loosely track ancestral diet type. Indicative "
                "only \u2014 not a clinical result.",
        "markers": {
            "rs4988235": "MCM6/LCT \u2014 adult lactase persistence ('farmer' "
                         "dairy adaptation)",
            "rs174547": "FADS1 \u2014 efficiency of making omega-3 from plants",
            "rs1801133": "MTHFR \u2014 folate handling on a plant-heavy diet",
            "rs1061325": "CLTCL1 \u2014 a diet-adaptation marker",
        },
    },
}

# Build gene -> (article_label, topic) index. First feature wins.
FEATURE_OF_GENE = {}
FEATURE_BASELINE = {}
for _label, _topic, _genes, _note in FEATURES:
    for _g in _genes:
        gu = _g.upper()
        FEATURE_OF_GENE.setdefault(gu, (_label, _topic))
        FEATURE_BASELINE.setdefault(gu, _note)

ARCHETYPE_MARKERS = {}            # rsid -> (archetype_label, descriptor)
for _alabel, _info in ARCHETYPES.items():
    for _rsid, _desc in _info["markers"].items():
        ARCHETYPE_MARKERS[_rsid] = (_alabel, _desc)


def feature_for(gene):
    """Return (article_label, topic_key) for a gene, or None."""
    return FEATURE_OF_GENE.get((gene or "").strip().upper())


def feature_baseline(gene):
    return FEATURE_BASELINE.get((gene or "").strip().upper(), "")


def is_feature_gene(gene):
    return (gene or "").strip().upper() in FEATURE_OF_GENE


# Register feature genes into gene_labels so they're treated as "known"
# (kept as their own articles) and so notes.py can fall back to the feature
# baseline. Done at import time.
try:
    from . import gene_labels as _gl
except ImportError:
    import gene_labels as _gl
for _label, _topic, _genes, _note in FEATURES:
    _short = _label.split(" (")[0]
    for _g in _genes:
        _gl.GENE_FUNCTION[_g.upper()] = (_short, _note)
