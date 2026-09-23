#!/usr/bin/env python3
"""ACMG-2015 automated rule engine, two arms.
Baseline arm: faithful re-implementation of InterVar's published automated tier
(Wang et al., J Mol Diagn 2020, PMC7257575): population-frequency criteria,
LoF criterion, in-silico concordance, same-codon ClinVar evidence.
Full arm: baseline + gnomAD constraint-informed criteria (LOEUF-modulated PVS1,
mis_z PP2/BP1). Every fired code carries its evidence datum (G4 traceability).
Combining rules: Richards et al. 2015 (PMID 25741868)."""
import csv, json, os, re

BASE = "/home/sandbox/rare-disease-variant-pipeline"
LOF_CSQ = {"transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
           "stop_gained", "frameshift_variant", "start_lost"}
LOF_GENES_CACHE = os.path.join(BASE, "data/cache/lof_genes.json")

def build_lof_genes():
    """Genes with >=3 distinct P/LP LoF variants in ClinVar (LoF-mechanism proxy,
    as InterVar uses a LoF-disease-gene list)."""
    if os.path.exists(LOF_GENES_CACHE):
        return set(json.load(open(LOF_GENES_CACHE)))
    pat = re.compile(r"(fs|Ter|\*)")
    counts = {}
    with open(os.path.join(BASE, "data/clinvar_grch38_small.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["ClinicalSignificance"] in ("Pathogenic", "Likely pathogenic",
                                               "Pathogenic/Likely pathogenic"):
                if pat.search(row["Name"]):
                    counts.setdefault(row["GeneSymbol"], set()).add(row["Name"])
    genes = {g for g, s in counts.items() if len(s) >= 3}
    json.dump(sorted(genes), open(LOF_GENES_CACHE, "w"))
    return genes

def classify(feat, lof_genes, arm="baseline"):
    """feat: merged feature dict. Returns (label, codes) where codes maps
    ACMG code -> evidence string. Label in P/LP/VUS/LB/B."""
    codes = {}
    csq = set(feat.get("consequence") or [])
    if not csq and feat.get("most_severe"):
        csq = {feat["most_severe"]}
    gene = feat.get("gene", "")
    af = feat.get("af")
    af = af if af is not None else 0.0
    sift = feat.get("sift", "")
    pp2_pred = feat.get("polyphen", "")
    is_missense = "missense_variant" in csq
    is_lof = bool(csq & LOF_CSQ)

    # --- population frequency (both arms; InterVar thresholds on gnomAD max AF) ---
    if af >= 0.05:
        codes["BA1"] = f"gnomAD max AF {af:.4f} >= 0.05"
    elif af >= 0.01:
        codes["BS1"] = f"gnomAD max AF {af:.4f} >= 0.01"
    elif af < 5e-5:
        codes["PM2"] = f"gnomAD max AF {af:.6f} < 5e-5 (absent/extremely rare)"

    # --- LoF ---
    if is_lof:
        if arm == "full":
            loeuf = feat.get("oe_lof_upper")
            if gene in lof_genes and loeuf is not None and loeuf < 0.35:
                codes["PVS1"] = f"LoF consequence in LoF-mechanism gene {gene} (LOEUF upper {loeuf:.3f}<0.35)"
            elif gene in lof_genes:
                codes["PVS1_Moderate"] = f"LoF in LoF-mechanism gene {gene} but LOEUF upper {loeuf if loeuf is None else round(loeuf,3)} >=0.35"
        else:
            if gene in lof_genes:
                codes["PVS1"] = f"LoF consequence in LoF-mechanism gene {gene} (ClinVar LoF proxy)"

    # --- same-codon evidence (both arms, other variants only) ---
    if feat.get("ps1"):
        codes["PS1"] = f"same amino-acid change previously P/LP in ClinVar (other record)"
    elif feat.get("pm5") and is_missense:
        codes["PM5"] = f"different P/LP missense change at same codon in ClinVar (other record)"

    # --- in-silico concordance (missense only) ---
    if is_missense:
        s_benign = sift == "tolerated" or sift.startswith("tolerated")
        s_del = sift.startswith("deleterious")
        p_benign = pp2_pred == "benign"
        p_del = pp2_pred in ("probably_damaging", "possibly_damaging")
        if sift and pp2_pred:
            if s_del and p_del:
                codes["PP3"] = f"SIFT={sift}, PolyPhen={pp2_pred} (concordant deleterious)"
            elif s_benign and p_benign:
                codes["BP4"] = f"SIFT={sift}, PolyPhen={pp2_pred} (concordant benign)"
        # BP1 / PP2
        if arm == "full":
            mis_z = feat.get("mis_z")
            if mis_z is not None and mis_z >= 3.09:
                codes["PP2"] = f"missense in missense-constrained gene {gene} (mis_z={mis_z:.2f}>=3.09)"
            elif mis_z is not None and mis_z < 0:
                codes["BP1"] = f"missense in missense-tolerant gene {gene} (mis_z={mis_z:.2f}<0)"
        else:
            if gene in lof_genes:
                codes["BP1"] = f"missense variant in gene {gene} with truncating mechanism"

    label = combine(codes)
    return label, codes

VS = {"PVS1"}
ST = {"PS1"}
MO = {"PM2", "PM5", "PVS1_Moderate"}
SU = {"PP2", "PP3"}
BS = {"BA1"}          # stand-alone
BST = {"BS1"}
BSU = {"BP1", "BP4"}

def combine(codes):
    n = lambda s: sum(1 for c in codes if c in s)
    vs, st, mo, su = n(VS), n(ST), n(MO), n(SU)
    ba, bs, bu = n(BS), n(BST), n(BSU)
    # benign first (stand-alone)
    if ba >= 1:
        return "B"
    if bs >= 2:
        return "B"
    if bs == 1 and bu >= 1:
        return "LB"
    if bu >= 2 and vs == 0 and st == 0:
        return "LB"
    # pathogenic
    if (vs >= 1 and (st >= 1 or mo >= 2 or (mo == 1 and su >= 2) or su >= 3)) or st >= 2 \
       or (st == 1 and (mo >= 3 or (mo == 2 and su >= 2) or (mo == 1 and su >= 4))):
        return "P"
    if (vs == 1 and (mo >= 1 or su >= 2)) or (st == 1 and (mo >= 1 or su >= 2)) \
       or mo >= 3 or (mo == 2 and su >= 2) or (mo == 1 and su >= 4):
        return "LP"
    return "VUS"

SCORE = {"P": 5, "LP": 4, "VUS": 3, "LB": 2, "B": 1}

if __name__ == "__main__":
    lg = build_lof_genes()
    print("LoF-mechanism genes:", len(lg), "e.g.", sorted(lg)[:5])
    # sanity: TP53 R248Q-like feature
    f = {"consequence": ["missense_variant"], "gene": "TP53", "af": 1.97e-5,
         "sift": "deleterious", "polyphen": "probably_damaging", "ps1": False,
         "pm5": True, "mis_z": 3.1, "oe_lof_upper": 0.9}
    for arm in ("baseline", "full"):
        print(arm, classify(f, lg, arm))
    # common benign
    f2 = {"consequence": ["missense_variant"], "gene": "TTN", "af": 0.2,
          "sift": "tolerated", "polyphen": "benign"}
    print("baseline", classify(f2, lg, "baseline"))
    # LoF
    f3 = {"consequence": ["frameshift_variant"], "gene": "BRCA2", "af": 0.0,
          "oe_lof_upper": 0.82}
    print("baseline", classify(f3, lg, "baseline"))
    print("full", classify(f3, lg, "full"))
