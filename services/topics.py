"""
Topic taxonomy for the GenoDynLabs condensed report.

This is OUR OWN categorization of variants into physiological systems. It is
NOT a copy of any third-party report's article structure. Two inputs drive it:

  1. A named GWAS trait (from trait_df.csv), matched by keyword -> topic.
  2. A gene symbol (from SNPedia), matched against a gene -> topic map.

Trait match wins; gene map is the fallback; anything left over goes to "Other".

The ordering of TOPICS controls the order sections appear in the report, and
each topic carries a pastel banner color (mirrors the look of the reference
layout without reusing its content).
"""

# (topic key, display label, banner background, banner text color)
TOPICS = [
    ("nutrients",   "Nutrients \u2013 Vitamins, Minerals & Diet", "#e6efd8", "#5c6b3f"),
    ("methylation", "Methylation & Homocysteine",                 "#dcefe4", "#2f6b4a"),
    ("diet",        "Diet, Food & Intolerances",                  "#f3e8d2", "#7a5c1e"),
    ("heart",       "Heart & Vascular Health",                    "#fff3d6", "#7a5c1e"),
    ("metabolic",   "Metabolic Health",                           "#d9e8f5", "#2f5169"),
    ("brain",       "Brain, Mood & Sleep",                        "#fbe0ea", "#7a3350"),
    ("immune",      "Immune, Allergy & Autoimmune",               "#e6e0f2", "#4a3d6b"),
    ("cancer",      "Cancer Risk",                                "#f6ddd6", "#7a3a28"),
    ("detox",       "Detoxification & Drug Metabolism",           "#d6ece8", "#235049"),
    ("hormones",    "Hormones & Fertility",                       "#efe7d6", "#6b5320"),
    ("musculo",     "Bone, Joint & Muscle",                       "#e0eaf2", "#33536b"),
    ("conditions",  "Inherited Conditions",                       "#f0e2e2", "#6b3a3a"),
    ("traits",      "Traits & Pharmacogenomics",                  "#f2ead6", "#6b5a20"),
    ("longevity",   "Longevity & Aging",                          "#ece0f2", "#523d6b"),
    ("other",       "Other Findings",                             "#ededed", "#444444"),
]

TOPIC_ORDER = [t[0] for t in TOPICS]
TOPIC_META = {t[0]: {"label": t[1], "bg": t[2], "fg": t[3]} for t in TOPICS}

# Keyword -> topic for matching named GWAS traits. Checked in order; first hit wins.
_TRAIT_RULES = [
    ("cancer",    ["carcinoma", "cancer", "melanoma", "lymphoma", "leukemia",
                   "neoplasm", "glioma", "myeloma", "sarcoma", "tumor", "adenoma"]),
    ("heart",     ["coronary", "myocardial", "atrial fibrillation", "heart",
                   "cardiac", "hypertension", "blood pressure", "stroke",
                   "thromboembolism", "aneurysm", "ldl", "hdl", "cholesterol",
                   "triglyceride", "lipoprotein", "vascular", "aortic"]),
    ("metabolic", ["diabetes", "glucose", "insulin", "obesity", "body mass",
                   "bmi", "metabolic", "fatty liver", "nafld", "waist",
                   "adiponectin", "uric acid", "gout"]),
    ("immune",    ["lupus", "rheumatoid", "crohn", "ulcerative", "inflammatory bowel",
                   "multiple sclerosis", "psoriasis", "asthma", "celiac",
                   "type 1 diabetes", "allergic", "allergy", "eczema",
                   "graves", "hashimoto", "sjogren", "autoimmune", "ankylosing",
                   "vitiligo", "ige", "atopic"]),
    ("brain",     ["schizophrenia", "depress", "bipolar", "alzheimer", "parkinson",
                   "anxiety", "cognitive", "intelligence", "neuroticism",
                   "migraine", "epilepsy", "sleep", "insomnia", "autism", "adhd",
                   "mood", "addiction", "alcohol depend", "smoking", "memory"]),
    ("musculo",   ["bone mineral", "osteoporosis", "osteoarthritis", "arthritis",
                   "height", "rheumatoid arthritis", "fracture", "hernia",
                   "tendin", "muscle", "back pain"]),
    ("hormones",  ["pcos", "polycystic", "endometri", "fibroid", "testosterone",
                   "estrogen", "menarche", "menopause", "fertility", "thyroid",
                   "prostate", "breast", "miscarriage", "preeclampsia"]),
    ("nutrients", ["vitamin", "folate", "homocysteine", "iron", "ferritin",
                   "selenium", "b12", "magnesium", "calcium", "zinc"]),
    ("longevity", ["longevity", "lifespan", "telomere", "aging", "age at death",
                   "centenarian", "frailty"]),
    ("detox",     ["drug", "warfarin", "metabolizer", "caffeine"]),
]

# Gene symbol -> topic, used when there is no named trait. Extend freely.
_GENE_TOPIC = {}
def _add(topic, *genes):
    for g in genes:
        _GENE_TOPIC[g.upper()] = topic

_add("nutrients", "MTHFR", "MTRR", "MTR", "FUT2", "TCN1", "TCN2", "FMO3", "BCMO1",
     "CYP2R1", "GC", "VDR", "DHFR", "SLC23A1", "BTD", "SLC30A8", "ALPL", "CBS",
     "MTHFD1", "SLC19A2", "PEMT", "CHKA", "BHMT", "SEP15", "SEPP1", "GPX1", "GPX4",
     "TRPM6", "CNNM2", "HFE", "HFE2", "SLC40A1", "TF", "TMPRSS6", "FADS1", "FADS2")
_add("heart", "PCSK9", "APOB", "LDLR", "LDLRAP1", "LPA", "CETP", "APOA5", "LPL",
     "F2", "F5", "F11", "FGA", "FGB", "FGG", "ITGB3", "NOS3", "CYP11B2", "AGTR1",
     "ACE", "SERPINE1", "SERPINA1", "PITX2", "ZFHX3", "ABCA1", "ABCG5", "ABCG8",
     "CRP", "GUCY1A3", "LIPC", "CDKN2B-AS1", "VWF", "MYBPC3", "MYH7", "TNNT2")
_add("metabolic", "TCF7L2", "MTNR1B", "PPARG", "PPARGC1A", "KCNJ11", "ABCC8",
     "GCK", "CDKAL1", "IRS1", "ENPP1", "GLP1R", "PCSK1", "SLC2A9", "ABCG2",
     "FTO", "MC4R", "LEPR", "ADIPOQ", "UCP1", "UCP2", "UCP3", "GHRL", "SH2B1",
     "PNPLA3", "TM6SF2", "MBOAT7", "GCKR", "SIRT1", "SIRT3", "SIRT6", "SCD",
     "SLC30A2", "NPY", "PPM1K", "HNF1A", "DHCR7")
_add("brain", "COMT", "BDNF", "OXTR", "OXT", "MAOA", "MAOB", "DRD1", "DRD2",
     "DRD3", "DRD4", "ANKK1", "HTR1A", "HTR1B", "HTR2A", "HTR2B", "HTR3A", "TPH2",
     "SLC6A4", "SLC6A3", "SLC6A2", "FKBP5", "CRHR1", "NR3C1", "NR3C2", "GAD1",
     "GAD2", "GABRA1", "GABRA2", "GABRA6", "CLOCK", "PER1", "PER2", "PER3",
     "CRY1", "CRY2", "OPN4", "ADA", "MEIS1", "BTBD9", "GSK3B", "ANK3", "CNR1",
     "CNR2", "FAAH", "AKT1", "SNAP25", "KIAA0319", "DCDC2", "SNCA", "LRRK2",
     "APOE", "TREM2", "CLU", "PICALM", "MAPT", "GADL1", "ACCN1", "ACCN2")
_add("immune", "TNF", "IL6", "IL6R", "IL1B", "IL1A", "IL10", "IL13", "IL17A",
     "IL17F", "IL2", "IL2RA", "IL12B", "IL7R", "TNFRSF1A", "TNFRSF1B", "TNFAIP3",
     "PTPN22", "CTLA4", "STAT4", "IRF5", "IRF8", "HLA-DRB1", "HLA-DQB1", "HLA-DQA1",
     "HLA-C", "HLA-B", "HLA-A", "NOD2", "ATG16L1", "IRGM", "CCR5", "MBL2", "OAS1",
     "TLR1", "TLR2", "TLR3", "TLR4", "TLR7", "NLRP3", "CIAS1", "MEFV", "FCGR2A",
     "AIRE", "BLK", "IFIH1", "CD6", "CD58", "ICAM1", "INFG", "IFNG")
_add("cancer", "BRCA1", "BRCA2", "MSH3", "HOXB13", "TP53", "CASC8", "CASC17",
     "GATA3", "AHR", "SMAD7", "HIF1A", "NQO1", "XPC", "XPG")
_add("detox", "CYP1A1", "CYP1A2", "CYP2A6", "CYP2B6", "CYP2C8", "CYP2C9",
     "CYP2C19", "CYP2D6", "CYP3A4", "CYP3A5", "CYP2E1", "NAT1", "NAT2", "GSTM1",
     "GSTP1", "GSTA1", "GSTO1", "SOD1", "SOD2", "UGT1A1", "UGT1A6", "NFE2L2",
     "ABCB1", "BCHE", "PON1", "AS3MT", "SLCO1B1", "CES1", "OPRM1", "SLC22A1")
_add("hormones", "CYP17A1", "CYP19A1", "CYP1B1", "ESR1", "ESR2", "SHBG", "PGR",
     "TSHR", "PDE8B", "FOXE1", "TPO", "DIO1", "DIO2", "LHCGR", "DENND1A", "FSHB",
     "FSHR", "SULT2A1", "SULT1A1", "HSD17B1", "HSD11B1", "PRM1", "GPER1")
_add("musculo", "ACTN3", "MSTN", "AMPD1", "COL1A1", "COL5A1", "COL3A1", "COL2A1",
     "COL11A1", "TNFSF11", "OPG", "LRP5", "VDR", "TGFB1", "ESR2", "MMP1", "MMP3",
     "MMP13", "GDF5", "PPARD", "AGT", "AGTR2")
_add("traits", "ABCC11", "MC1R", "OR5A1", "OR7D4", "TAS2R38", "TAS2R16",
     "TAS1R2", "TAS1R3", "ALDH2", "ADH1B", "ADH1C", "LCT", "MCM6", "CYP1A2",
     "TYK2", "AMY1", "SLC22A4", "DEC2")
_add("longevity", "FOXO3A", "IGF1R", "KL", "TERT", "TERC", "NAF1", "OBFC1",
     "NAMPT", "CD38", "BST1", "AGER", "GLO1", "HSPA1L", "HSPA5", "TRAP1",
     "MTOR", "RPTOR", "DEPTOR", "BPIFB4")


# Curated genes that were missing from the maps above. Keeps everything in
# gene_labels.GENE_FUNCTION out of the "Other" bucket.
_add("nutrients", "FOLH1", "SLC19A1", "GGCX", "TFR2", "SELENOP", "APOA2")
_add("heart", "9p21", "VKORC1")
_add("detox", "TPMT", "DPYD")
_add("brain", "BHLHE41")
_add("immune", "IL23R")
_add("longevity", "FOXO3", "KLOTHO")

# Route every nutrient-panel gene to the Nutrients topic so vitamins/minerals
# always land in the Vitamins & Minerals section.
try:
    from . import gene_labels as _gl
except ImportError:
    import gene_labels as _gl
_add("nutrients", *_gl.NUTRIENT_OF_GENE.keys())


def classify_trait(trait):
    if not trait or str(trait).strip().lower() in ("", "nan", "none", "nr"):
        return None
    t = str(trait).lower()
    for topic, kws in _TRAIT_RULES:
        if any(k in t for k in kws):
            return topic
    return None


def classify_gene(gene):
    if not gene:
        return None
    g = str(gene).strip().upper()
    # SNPedia sometimes lists several genes joined by commas/spaces; try each.
    for token in g.replace(";", ",").replace(" ", ",").split(","):
        token = token.strip()
        if token in _GENE_TOPIC:
            return _GENE_TOPIC[token]
    return None


def classify(gene=None, trait=None):
    """Return a topic key. Trait keyword match wins, then gene map, then 'other'."""
    return classify_trait(trait) or classify_gene(gene) or "other"
