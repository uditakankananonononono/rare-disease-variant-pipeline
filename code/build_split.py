#!/usr/bin/env python3
"""Build variant-level tune/held-out split from ClinVar GRCh38 small-variant subset.
G5: split by variant (md5 of VariationID), deterministic. No ClinVar label leaks
into features downstream; labels used only as answer key."""
import csv, hashlib, random, collections

CLINVAR = "/home/sandbox/rare-disease-variant-pipeline/data/clinvar_grch38_small.tsv"
G2P = "/home/sandbox/rare-disease-variant-pipeline/data/genes_to_phenotype.txt"
OUT = "/home/sandbox/rare-disease-variant-pipeline/repo/data/variant_eval_split.tsv"

GOOD_REVIEW = {"reviewed by expert panel", "practice guideline",
               "criteria provided, multiple submitters, no conflicts"}
POS = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
NEG = {"Benign", "Likely benign", "Benign/Likely benign"}

mendelian = set()
with open(G2P) as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        mendelian.add(row["gene_symbol"])

seen = set()
cands = []
with open(CLINVAR) as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        sig = row["ClinicalSignificance"]
        if sig not in POS and sig not in NEG:
            continue
        if row["ReviewStatus"] not in GOOD_REVIEW:
            continue
        if row["GeneSymbol"] not in mendelian:
            continue
        if not row["PositionVCF"] or row["PositionVCF"] == "na":
            continue
        key = (row["Chromosome"], row["PositionVCF"], row["ReferenceAlleleVCF"], row["AlternateAlleleVCF"])
        if key in seen or "na" in key:
            continue
        seen.add(key)
        label = 1 if sig in POS else 0
        h = int(hashlib.md5(row["VariationID"].encode()).hexdigest(), 16)
        cands.append({**row, "label": label, "split": "heldout" if h % 10 >= 8 else "tune"})

random.seed(42)
pos = [c for c in cands if c["label"] == 1]
neg = [c for c in cands if c["label"] == 0]
random.shuffle(pos); random.shuffle(neg)
N = 3000
sel = pos[:N] + neg[:N]
random.shuffle(sel)

with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["variation_id", "chrom", "pos", "ref", "alt", "gene", "sig", "label",
                "review_status", "name", "phenotype_ids", "split"])
    for c in sel:
        w.writerow([c["VariationID"], c["Chromosome"], c["PositionVCF"],
                    c["ReferenceAlleleVCF"], c["AlternateAlleleVCF"], c["GeneSymbol"],
                    c["ClinicalSignificance"], c["label"], c["ReviewStatus"],
                    c["Name"], c["PhenotypeIDS"], c["split"]])

print("candidate pool:", len(cands), "pos:", len(pos), "neg:", len(neg))
print("selected:", len(sel))
print(collections.Counter((c["split"], c["label"]) for c in sel))
