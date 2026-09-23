#!/usr/bin/env python3
"""Annotate the 6,000-variant eval split: fetch features (disk-cached), classify
both arms, write results TSV. Resume-safe."""
import csv, json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from features import vep_batch, gnomad_batch, constraint_batch, build_codon_index, same_codon, vid
from acmg import classify, build_lof_genes, SCORE

BASE = "/home/sandbox/rare-disease-variant-pipeline"
SPLIT = os.path.join(BASE, "repo/data/variant_eval_split.tsv")
OUT = os.path.join(BASE, "results_variant_level.tsv")
LOG = os.path.join(BASE, "annotate.log")

def log(msg):
    with open(LOG, "a") as fh:
        fh.write(time.strftime("%H:%M:%S ") + msg + "\n")

rows = list(csv.DictReader(open(SPLIT), delimiter="\t"))
log(f"loaded {len(rows)} split variants")
variants = [(r["chrom"], r["pos"], r["ref"], r["alt"]) for r in rows]

t0 = time.time()
vep = vep_batch(variants)
log(f"vep done {time.time()-t0:.0f}s")
gn = gnomad_batch(variants)
log(f"gnomad done {time.time()-t0:.0f}s")
genes = sorted({v.get("gene", "") for v in vep.values() if v.get("gene")})
con = constraint_batch(genes)
log(f"constraint done {len(genes)} genes {time.time()-t0:.0f}s")
ci = build_codon_index()
lg = build_lof_genes()
log("indexes built")

with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["variation_id", "vid", "gene", "label", "split", "sig",
                "most_severe", "af", "sift", "polyphen", "ps1", "pm5",
                "oe_lof_upper", "mis_z", "baseline_label", "baseline_codes",
                "full_label", "full_codes"])
    for r in rows:
        k = vid(r["chrom"], r["pos"], r["ref"], r["alt"])
        v = vep.get(k, {})
        g = gn.get(k, {"af": None})
        cg = con.get(v.get("gene", ""), {})
        ps1, pm5 = same_codon(ci, v.get("gene", ""), v.get("protein_start"),
                              v.get("amino_acids", ""), k)
        feat = {**v, **g, **cg, "ps1": ps1, "pm5": pm5}
        bl, bc = classify(feat, lg, "baseline")
        fl, fc = classify(feat, lg, "full")
        w.writerow([r["variation_id"], k, v.get("gene", ""), r["label"], r["split"],
                    r["sig"], v.get("most_severe", ""), g.get("af"),
                    v.get("sift", ""), v.get("polyphen", ""), ps1, pm5,
                    cg.get("oe_lof_upper"), cg.get("mis_z"), bl,
                    json.dumps(bc), fl, json.dumps(fc)])
log(f"DONE {time.time()-t0:.0f}s -> {OUT}")
