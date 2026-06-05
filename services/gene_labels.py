"""
gene_labels.py  --  Curated gene -> biological-function labels and baseline notes.

WHY THIS EXISTS
The report used to title every article with the raw gene symbol (e.g. "FUT2").
That is opaque to a user and, in context, wrong: in the Nutrients section a FUT2
variant is really about vitamin B12, so the article should read "Vitamin B12
(FUT2)". This module is the single source of truth for:

  1. label(gene)         -> a human, function-first article title
  2. baseline_note(gene) -> a one-line "what this gene does" sentence used as a
                            LAST-RESORT note so no gene pair is ever blank.

The map is hand-curated from well-established gene/function relationships. It is
intentionally conservative: entries describe the gene's primary, textbook role,
not speculative associations. Anything not in the map falls back gracefully to
the strongest GWAS trait (handled in notes.py) and finally to the gene symbol.

Add freely; keys are upper-case gene symbols.
"""

# gene -> (short function label, one-line baseline description)
# The label is what shows in the article header *before* the "(GENE)" suffix.
GENE_FUNCTION = {
    # ---- B-vitamins, folate, methylation -------------------------------
    "MTHFR":  ("Folate Metabolism", "Converts folate to its active form; variants slow folic-acid processing and can raise homocysteine."),
    "MTRR":   ("Folate / B12 Recycling", "Regenerates the B12 cofactor that drives the methionine cycle."),
    "MTR":    ("B12-Dependent Methylation", "Uses vitamin B12 to recycle homocysteine back to methionine."),
    "MTHFD1": ("Folate One-Carbon Pool", "Channels folate through the one-carbon pathway used for DNA synthesis."),
    "FUT2":   ("Vitamin B12 Status", "Secretor gene; strongly influences serum vitamin B12 and gut microbiome composition."),
    "TCN1":   ("B12 Transport", "Encodes a B12-binding transport protein (haptocorrin)."),
    "TCN2":   ("B12 Transport", "Encodes transcobalamin, the protein that delivers B12 to cells."),
    "CBS":    ("Homocysteine Clearance", "Clears homocysteine via the transsulfuration pathway."),
    "BHMT":   ("Choline / Homocysteine", "Uses choline-derived betaine to remethylate homocysteine."),
    "PEMT":   ("Choline Synthesis", "Synthesizes phosphatidylcholine; variants raise dietary choline needs."),
    "FOLH1":  ("Dietary Folate Uptake", "Cleaves dietary folate so it can be absorbed in the gut."),
    "DHFR":   ("Folate Activation", "Reduces folic acid to its usable tetrahydrofolate form."),
    "SLC19A1":("Folate Transport", "Moves folate into cells."),

    # ---- Fat-soluble & other vitamins ----------------------------------
    "GC":     ("Vitamin D Transport", "Vitamin-D binding protein; influences circulating vitamin D."),
    "VDR":    ("Vitamin D Response", "Vitamin D receptor; mediates vitamin D's effects on bone and immunity."),
    "CYP2R1": ("Vitamin D Activation", "Hydroxylates vitamin D to its circulating form."),
    "BCMO1":  ("Vitamin A from Beta-Carotene", "Splits dietary beta-carotene into vitamin A; variants reduce conversion."),
    "GGCX":   ("Vitamin K Use", "Uses vitamin K to activate clotting and bone proteins."),

    # ---- Minerals ------------------------------------------------------
    "HFE":    ("Iron Overload (Hemochromatosis)", "Regulates iron absorption; variants can cause iron overload."),
    "TMPRSS6":("Iron Regulation", "Controls hepcidin, the master switch for iron absorption."),
    "TF":     ("Iron Transport", "Transferrin; carries iron in the blood."),
    "SLC40A1":("Iron Export", "Ferroportin; exports iron from cells."),
    "TFR2":   ("Iron Sensing", "Senses body iron stores to tune absorption."),
    "SLC30A8":("Zinc & Insulin", "Zinc transporter in insulin-secreting beta cells."),
    "TRPM6":  ("Magnesium Handling", "Reabsorbs magnesium in the kidney and gut."),
    "SELENOP":("Selenium Transport", "Carries selenium to tissues."),
    "GPX1":   ("Selenium Antioxidant", "Selenium-dependent antioxidant enzyme."),

    # ---- Fatty acids ---------------------------------------------------
    "FADS1":  ("Omega-3/6 Processing", "Desaturase that builds long-chain omega-3 and omega-6 fats."),
    "FADS2":  ("Omega-3/6 Processing", "Desaturase that builds long-chain omega-3 and omega-6 fats."),
    "APOA2":  ("Saturated-Fat Response", "Variants link saturated-fat intake to weight gain."),

    # ---- Lipids / cardiovascular --------------------------------------
    "APOE":   ("Alzheimer's & Cholesterol Risk", "Carries cholesterol in blood and brain; the strongest common genetic factor for late-onset Alzheimer's."),
    "PCSK9":  ("LDL Cholesterol", "Controls how many LDL receptors clear cholesterol from blood."),
    "LDLR":   ("LDL Cholesterol Clearance", "The receptor that removes LDL cholesterol from blood."),
    "APOB":   ("LDL Particle Levels", "Main protein of LDL ('bad') cholesterol particles."),
    "LPA":    ("Lipoprotein(a)", "Sets lifelong Lp(a), an independent cardiovascular risk factor."),
    "CETP":   ("HDL Cholesterol", "Transfers cholesterol between HDL and other particles."),
    "LPL":    ("Triglyceride Clearance", "Breaks down triglycerides in the bloodstream."),
    "F5":     ("Clotting (Factor V Leiden)", "Factor V; the Leiden variant raises clot risk."),
    "F2":     ("Clotting (Prothrombin)", "Prothrombin; variants raise venous clot risk."),
    "NOS3":   ("Blood-Vessel Tone", "Makes nitric oxide that relaxes blood vessels."),
    "ACE":    ("Blood Pressure", "Angiotensin-converting enzyme; regulates blood pressure."),
    "AGT":    ("Blood Pressure", "Angiotensinogen; precursor in blood-pressure control."),
    "9p21":   ("Coronary Artery Disease", "The most replicated common-variant locus for heart disease."),
    "CDKN2B-AS1": ("Coronary Artery Disease", "9p21 locus; the most replicated heart-disease region."),

    # ---- Glucose / metabolic ------------------------------------------
    "TCF7L2": ("Type 2 Diabetes Risk", "Strongest common-variant risk gene for type 2 diabetes."),
    "FTO":    ("Body Weight / Obesity", "Best-known common variant influencing BMI and appetite."),
    "MC4R":   ("Appetite & Obesity", "Hypothalamic switch controlling appetite and body weight."),
    "PPARG":  ("Insulin Sensitivity", "Master regulator of fat cells and insulin sensitivity."),
    "MTNR1B": ("Fasting Glucose", "Melatonin receptor linking sleep timing to glucose control."),
    "SLC2A9": ("Uric Acid / Gout", "Kidney urate transporter; sets uric-acid levels."),
    "ABCG2":  ("Uric Acid / Gout", "Excretes uric acid; variants raise gout risk."),
    "PNPLA3": ("Fatty Liver", "Strongest common variant for non-alcoholic fatty liver."),

    # ---- Lactose, alcohol, caffeine, taste ----------------------------
    "MCM6":   ("Lactose Tolerance", "Regulatory region controlling adult lactase persistence."),
    "LCT":    ("Lactose Tolerance", "Lactase enzyme; determines adult milk digestion."),
    "ALDH2":  ("Alcohol Flush", "Clears acetaldehyde; the deficient variant causes the alcohol flush."),
    "ADH1B":  ("Alcohol Metabolism", "Speed of converting alcohol to acetaldehyde."),
    "CYP1A2": ("Caffeine Metabolism", "Main enzyme that clears caffeine; sets fast vs slow status."),
    "TAS2R38":("Bitter Taste", "Determines whether you taste certain bitter compounds (PTC/PROP)."),
    "TAS1R2": ("Sweet Preference", "Sweet-taste receptor."),
    "ABCC11": ("Earwax & Body Odor", "Determines wet vs dry earwax and underarm odor."),

    # ---- Detox / pharmacogenomics -------------------------------------
    "CYP2C19":("Drug Metabolism (CYP2C19)", "Activates/clears clopidogrel, PPIs, and many antidepressants."),
    "CYP2D6": ("Drug Metabolism (CYP2D6)", "Clears ~25% of common drugs including codeine and many antidepressants."),
    "CYP2C9": ("Warfarin Sensitivity", "Clears warfarin and NSAIDs; sets dosing."),
    "VKORC1": ("Warfarin Sensitivity", "Warfarin's target; the main dosing determinant."),
    "SLCO1B1":("Statin Side-Effects", "Liver transporter; variants raise statin muscle-pain risk."),
    "TPMT":   ("Thiopurine Drugs", "Clears thiopurine chemo/immune drugs; deficiency causes toxicity."),
    "DPYD":   ("Fluoropyrimidine Chemo", "Clears 5-FU/capecitabine; deficiency causes severe toxicity."),
    "UGT1A1": ("Bilirubin (Gilbert's)", "Clears bilirubin and irinotecan; variants cause Gilbert's syndrome."),
    "NAT2":   ("Acetylator Status", "Fast vs slow acetylation of drugs and dietary amines."),
    "GSTM1":  ("Toxin Detox (GST)", "Glutathione transferase; detoxifies environmental chemicals."),
    "COMT":   ("Dopamine Breakdown", "Clears dopamine in the prefrontal cortex; affects stress and focus."),

    # ---- Brain / mood / sleep -----------------------------------------
    "BDNF":   ("Memory & Mood (BDNF)", "Growth factor for neurons; tied to memory and mood."),
    "OXTR":   ("Social Bonding", "Oxytocin receptor; linked to empathy and social behavior."),
    "MAOA":   ("Mood & Aggression", "Breaks down serotonin, dopamine, and noradrenaline."),
    "SLC6A4": ("Serotonin Transport", "The serotonin transporter targeted by SSRIs."),
    "CLOCK":  ("Circadian Rhythm", "Core clock gene influencing chronotype and sleep timing."),
    "PER3":   ("Chronotype / Sleep", "Clock gene linked to morning vs evening preference."),
    "DEC2":   ("Short Sleep", "Rare variants allow healthy short sleep (BHLHE41)."),
    "BHLHE41":("Short Sleep", "Rare variants allow healthy short sleep (DEC2)."),

    # ---- Immune / autoimmune ------------------------------------------
    "HLA-DRB1": ("Autoimmune Risk (HLA)", "Immune presentation gene; central to autoimmune risk."),
    "HLA-DQB1": ("Celiac / Autoimmune (HLA)", "Immune gene; determines celiac and other autoimmune risk."),
    "PTPN22": ("Autoimmune Risk", "Immune-cell brake; broad autoimmune risk variant."),
    "NOD2":   ("Crohn's Disease", "Bacterial sensor; strongest common variant for Crohn's."),
    "IL23R":  ("Inflammatory Bowel Disease", "Inflammation receptor tied to IBD risk and protection."),
    "CCR5":   ("HIV Resistance", "HIV co-receptor; the delta-32 deletion confers resistance."),

    # ---- Cancer --------------------------------------------------------
    "BRCA1":  ("Hereditary Breast/Ovarian Cancer", "DNA-repair gene; pathogenic variants sharply raise breast/ovarian cancer risk."),
    "BRCA2":  ("Hereditary Breast/Ovarian Cancer", "DNA-repair gene; pathogenic variants raise breast/ovarian/prostate cancer risk."),
    "TP53":   ("Tumor Suppressor (TP53)", "Guardian of the genome; central tumor suppressor."),
    "MC1R":   ("Red Hair & Melanoma Risk", "Pigment gene; variants give red hair and raise melanoma risk."),

    # ---- Bone / muscle / performance ----------------------------------
    "ACTN3":  ("Muscle Fiber Type", "The 'sprinter gene'; affects fast-twitch muscle and power."),
    "COL1A1": ("Bone & Connective Tissue", "Builds type-I collagen for bone and tendon."),
    "COL5A1": ("Tendon & Ligament", "Type-V collagen; tied to tendon flexibility and injury risk."),

    # ---- Longevity -----------------------------------------------------
    "FOXO3":  ("Longevity", "One of the few replicated human longevity genes."),
    "FOXO3A": ("Longevity", "One of the few replicated human longevity genes."),
    "KLOTHO": ("Aging / Cognition", "Anti-aging hormone gene tied to cognition and lifespan."),
    "KL":     ("Aging / Cognition", "Klotho; anti-aging hormone tied to cognition and lifespan."),
}

# Some SNPedia gene fields use aliases; normalize a few common ones.
_ALIASES = {
    "FOXO3A": "FOXO3",
    "KLOTHO": "KL",
}


def _canon(gene):
    g = (gene or "").strip().upper()
    return _ALIASES.get(g, g)


def label(gene, fallback_trait=""):
    """Function-first article label.

    'FUT2'                      -> 'Vitamin B12 Status (FUT2)'
    'SOME_GENE' + trait 'Gout'  -> 'Gout (SOME_GENE)'
    'SOME_GENE' (nothing known) -> 'SOME_GENE'
    """
    g = _canon(gene)
    if g in GENE_FUNCTION:
        return f"{GENE_FUNCTION[g][0]} ({g})"
    if fallback_trait:
        ft = fallback_trait.strip()
        if ft and ft.lower() not in ("nan", "none", "nr"):
            ft = ft[:1].upper() + ft[1:]
            return f"{ft} ({g})" if g else ft
    return g or "Other variants"


def baseline_note(gene):
    """Last-resort, always-true sentence about what the gene does. Empty if the
    gene isn't in the curated map (notes.py then tries other sources)."""
    g = _canon(gene)
    info = GENE_FUNCTION.get(g)
    return info[1] if info else ""


def is_known(gene):
    return _canon(gene) in GENE_FUNCTION


# ====================================================================
# NUTRIENT PANEL
# --------------------------------------------------------------------
# The Vitamins & Minerals section is special: users expect to see a fixed
# panel of nutrients (not scattered gene articles). So we map established
# nutrient genes to a NUTRIENT, and the report groups all of a nutrient's
# genes under one article (e.g. "Vitamin B12" listing FUT2, TCN2, MUT...).
#
# Each entry: nutrient label -> (always-true baseline note, [canonical genes]).
# Genes not present in a given user's file simply don't appear. Genes are
# assigned to ONE primary nutrient (the most established association) so each
# lands in a single article; the Gene column preserves the real symbol.
#
# Ordered: fat-soluble + C, then B-complex, then macro/trace minerals.
# Grounded in established nutrigenomics (transporters, cofactor enzymes,
# absorption receptors). Add genes freely.
# ====================================================================
NUTRIENTS = [
    ("Vitamin A",
     "Influences how efficiently dietary beta-carotene is converted to active vitamin A and how retinol is carried in blood.",
     ["BCO1", "BCMO1", "BCO2", "RBP4", "TTR"]),
    ("Vitamin B1 (Thiamine)",
     "Affects thiamine transport into cells and its activation to the coenzyme form used in energy metabolism.",
     ["SLC19A2", "SLC19A3", "TPK1", "SLC25A19"]),
    ("Vitamin B2 (Riboflavin)",
     "Affects riboflavin transport and the FAD/FMN cofactors many enzymes (including MTHFR) depend on.",
     ["SLC52A1", "SLC52A2", "SLC52A3"]),
    ("Vitamin B3 (Niacin)",
     "Part of the pathway that makes NAD (the active form of niacin) from tryptophan and dietary niacin.",
     ["HAAO", "QPRT", "NAPRT", "NNMT", "NMNAT1", "NMNAT2", "NMNAT3", "NADK"]),
    ("Vitamin B5 (Pantothenic Acid)",
     "Affects pantothenate handling and coenzyme A synthesis used throughout fat and energy metabolism.",
     ["PANK1", "PANK2", "PANK3", "PANK4", "SLC5A6"]),
    ("Vitamin B6",
     "Affects activation and levels of vitamin B6 (pyridoxal-5'-phosphate), a cofactor for ~150 enzymes.",
     ["NBPF3", "ALPL", "ALDH7A1", "PNPO"]),
    ("Vitamin B7 (Biotin)",
     "Affects recycling and use of biotin; deficiency in these enzymes impairs the body's biotin economy.",
     ["BTD", "HLCS"]),
    ("Vitamin B9 (Folate)",
     "Affects folate transport and the folate pool; the MTHFR/folate-cycle enzymes are shown under Methylation & Folate Cycle.",
     ["DHFR", "SLC19A1", "FOLR1", "FOLH1", "GGH", "TYMS", "FPGS"]),
    ("Vitamin B12",
     "Affects vitamin B12 absorption and transport; the B12-dependent recycling enzymes are shown under Methylation & Folate Cycle.",
     ["FUT2", "FUT6", "TCN1", "TCN2", "CUBN", "AMN", "GIF", "MMAA", "MMAB", "MUT"]),
    ("Vitamin C",
     "Affects how vitamin C (ascorbate) is taken into cells and the plasma vitamin C level the body maintains.",
     ["SLC23A1", "SLC23A2"]),
    ("Vitamin D",
     "Affects vitamin D synthesis, transport, activation, and breakdown; a major driver of circulating vitamin D levels.",
     ["GC", "VDR", "CYP2R1", "CYP27B1", "CYP24A1", "DHCR7", "NADSYN1"]),
    ("Vitamin E",
     "Affects how vitamin E (tocopherol) is carried and distributed in blood.",
     ["TTPA", "APOA5", "SCARB1", "CD36", "ZNF259", "BUD13"]),
    ("Vitamin K",
     "Affects vitamin K recycling and the activation of clotting and bone proteins (and warfarin response).",
     ["GGCX", "VKORC1", "CYP4F2"]),
    ("Calcium",
     "Affects calcium sensing, absorption, and regulation in blood and bone.",
     ["CASR", "GCKR", "WDR81", "DGKD", "CARS"]),
    ("Magnesium",
     "Affects magnesium reabsorption and transport in the kidney and gut.",
     ["TRPM6", "TRPM7", "CNNM2", "SHROOM3", "ATP2B1", "SLC41A1", "MUC1"]),
    ("Sodium",
     "Affects sodium handling and salt sensitivity of blood pressure.",
     ["ADD1", "AGT", "CYP11B2", "SLC12A3", "GNB3", "NR3C2", "SCNN1B"]),
    ("Potassium",
     "Affects potassium handling in the kidney.",
     ["KCNJ1", "SLC12A1", "WNK1", "WNK4"]),
    ("Iron",
     "Affects iron absorption, transport, and storage; some variants raise the risk of iron overload.",
     ["HFE", "TMPRSS6", "TF", "TFR2", "TFRC", "SLC40A1", "HAMP", "BMP2", "FTL", "FTH1", "ARSB"]),
    ("Zinc",
     "Affects zinc transport and cellular zinc balance.",
     ["SLC30A8", "SLC30A3"]),
    ("Copper",
     "Affects copper transport and incorporation into ceruloplasmin; rare variants cause copper-handling disorders.",
     ["CP", "ATP7A", "ATP7B", "SLC31A1"]),
    ("Manganese",
     "Affects manganese transport; variants can shift manganese toward deficiency or overload.",
     ["SLC30A10", "SLC39A14", "SLC39A8"]),
    ("Iodine",
     "Affects iodine uptake and thyroid-hormone production that depends on it.",
     ["SLC5A5", "TPO", "TG", "DIO1", "DIO2"]),
    ("Selenium",
     "Affects selenium transport and the selenoprotein antioxidant enzymes that use it.",
     ["GPX1", "SELENOP", "SEPP1", "SELENBP1"]),
    ("Molybdenum",
     "Affects synthesis of the molybdenum cofactor required by sulfite oxidase and related enzymes.",
     ["MOCS1", "MOCS2", "GPHN", "SUOX"]),
]

# Chromium is the one nutrient with no gene to show: no human gene is robustly
# established to affect chromium status, so it stays an honest placeholder.
NUTRIENTS_LIMITED = {
    "Chromium":
        "No human gene is robustly established to affect chromium status, so there are no variants to report.",
}

# Build lookups and fold nutrient genes into GENE_FUNCTION (nutrient labels win).
NUTRIENT_ORDER = [n[0] for n in NUTRIENTS]
NUTRIENT_OF_GENE = {}
NUTRIENT_BASELINE = {}
for _label, _note, _genes in NUTRIENTS:
    NUTRIENT_BASELINE[_label] = _note
    for _g in _genes:
        gu = _g.upper()
        NUTRIENT_OF_GENE[gu] = _label
        # nutrient assignment is authoritative for label + baseline
        GENE_FUNCTION[gu] = (_label, _note)


def nutrient_of(gene):
    """Return the nutrient label for a gene, or None if it isn't a nutrient gene."""
    return NUTRIENT_OF_GENE.get(_canon(gene))


def is_nutrient_gene(gene):
    return _canon(gene) in NUTRIENT_OF_GENE


# The FULL fixed panel, in the exact order requested, so EVERY nutrient appears
# in every report — with the user's variants if any were found, else an honest
# "not assessed in your data" row. B3/Niacin and Chromium have no usable common
# variants; Vitamin C transporters are rarely on consumer chips.
PANEL_ORDER = [
    "Vitamin A", "Vitamin B1 (Thiamine)", "Vitamin B2 (Riboflavin)",
    "Vitamin B3 (Niacin)", "Vitamin B5 (Pantothenic Acid)", "Vitamin B6",
    "Vitamin B7 (Biotin)", "Vitamin B9 (Folate)", "Vitamin B12", "Vitamin C",
    "Vitamin D", "Vitamin E", "Vitamin K",
    "Calcium", "Magnesium", "Sodium", "Potassium", "Iron", "Zinc", "Copper",
    "Manganese", "Iodine", "Selenium", "Chromium", "Molybdenum",
]


def panel_empty_note(label):
    """Note shown when no variant for this nutrient was found in the user's data."""
    if label in NUTRIENTS_LIMITED:
        return NUTRIENTS_LIMITED[label]
    base = NUTRIENT_BASELINE.get(label, "")
    lead = "No variants affecting this nutrient were found in your uploaded data."
    return f"{lead} {base}".strip()
