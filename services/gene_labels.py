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
