#!/usr/bin/env python3
"""Parse phenopacket-store into a case index (data manifest artifact).
Extracts per-phenopacket: ids, disease, observed HPO terms, causal variant
(hg38 VCF record), allelic state, ACMG class. No outcome scoring here."""
import json, os, sys, csv

ROOT = sys.argv[1]
OUT = sys.argv[2]

rows = []
stats = {"total": 0, "with_vcf": 0, "hg38": 0, "small_variant": 0,
         "acmg_p_or_lp": 0, "hpo_ge3": 0, "evaluable": 0}

for cohort in sorted(os.listdir(ROOT)):
    pdir = os.path.join(ROOT, cohort, "phenopackets")
    if not os.path.isdir(pdir):
        continue
    for fn in sorted(os.listdir(pdir)):
        if not fn.endswith(".json"):
            continue
        stats["total"] += 1
        try:
            d = json.load(open(os.path.join(pdir, fn)))
        except Exception:
            continue
        pid = d.get("id", fn[:-5])
        hpo = [f["type"]["id"] for f in d.get("phenotypicFeatures", [])
               if not f.get("excluded", False) and "type" in f]
        diseases = d.get("diseases", [])
        dis_id = diseases[0]["term"]["id"] if diseases and "term" in diseases[0] else ""
        dis_lb = diseases[0]["term"].get("label", "") if diseases and "term" in diseases[0] else ""
        sex = d.get("subject", {}).get("sex", "")
        for it in d.get("interpretations", []):
            diag = it.get("diagnosis", {})
            for gi in diag.get("genomicInterpretations", []):
                vi = gi.get("variantInterpretation", {})
                acmg = vi.get("acmgPathogenicityClassification", "")
                vd = vi.get("variationDescriptor", {})
                gene = vd.get("geneContext", {}).get("valueId", "")
                allst = vd.get("allelicState", {}).get("label", "")
                hgvs_c = next((e["value"] for e in vd.get("expressions", [])
                               if e.get("syntax") == "hgvs.c"), "")
                hgvs_p = next((e["value"] for e in vd.get("expressions", [])
                               if e.get("syntax") == "hgvs.p"), "")
                vcf = vd.get("vcfRecord")
                if not vcf:
                    continue
                stats["with_vcf"] += 1
                asm = vcf.get("genomeAssembly", "")
                if asm not in ("hg38", "GRCh38"):
                    continue
                stats["hg38"] += 1
                chrom = (vcf.get("chrom") or "").replace("chr", "")
                pos = vcf.get("pos", 0)
                ref = vcf.get("ref", "") or ""
                alt = vcf.get("alt", "") or ""
                if chrom not in [str(c) for c in range(1, 23)] + ["X"]:
                    continue
                if max(len(ref), len(alt)) > 50:
                    continue
                stats["small_variant"] += 1
                if acmg in ("PATHOGENIC", "LIKELY_PATHOGENIC"):
                    stats["acmg_p_or_lp"] += 1
                if len(hpo) >= 3:
                    stats["hpo_ge3"] += 1
                eval_flag = (acmg in ("PATHOGENIC", "LIKELY_PATHOGENIC")
                             and len(hpo) >= 3 and gene)
                if eval_flag:
                    stats["evaluable"] += 1
                rows.append([pid, cohort, dis_id, dis_lb, sex, gene, chrom, pos,
                             ref, alt, allst, acmg, hgvs_c, hgvs_p,
                             len(hpo), ";".join(hpo), int(bool(eval_flag))])

with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["phenopacket_id", "cohort", "disease_id", "disease_label", "sex",
                "gene", "chrom", "pos", "ref", "alt", "allelic_state", "acmg",
                "hgvs_c", "hgvs_p", "n_hpo", "hpo_ids", "evaluable"])
    w.writerows(rows)

print(json.dumps(stats, indent=1))
print("rows:", len(rows))
