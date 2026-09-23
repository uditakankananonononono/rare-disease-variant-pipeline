#!/usr/bin/env python3
"""Select the 400 benchmark cases and write case shells (no background fetch;
run_cases.assemble_cases attaches HG002 background from the per-gene cache)."""
import csv, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_cases import build_real, build_sim

BASE = "/home/sandbox/rare-disease-variant-pipeline"
outdir = os.path.join(BASE, "data/cases")
os.makedirs(outdir, exist_ok=True)
real = build_real(200)
sim = build_sim(200)
with open(os.path.join(BASE, "data/cases_index.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["case_id", "arm", "causal_gene", "disease", "n_hpo", "panel_size",
                "causal_in_panel"])
    for c in real + sim:
        c["causal_in_panel"] = c["causal_gene"] in c["panel"]
        c["candidates"] = []
        json.dump(c, open(os.path.join(outdir, c["case_id"] + ".json"), "w"))
        w.writerow([c["case_id"], c["arm"], c["causal_gene"], c["disease"],
                    len(c["hpo"]), len(c["panel"]), c["causal_in_panel"]])
ri = sum(1 for c in real if c["causal_in_panel"])
si = sum(1 for c in sim if c["causal_in_panel"])
print(f"real: {len(real)} (causal-in-panel {ri}), sim: {len(sim)} (causal-in-panel {si})")
print("unique panel genes:", len({g for c in real + sim for g in c['panel']}))
