#!/usr/bin/env python3
"""G3/G4 evaluation from results/case_level_results.tsv.
G3: top-5 causal-variant recall (full arm) must exceed baseline arm by >=5pp
with bootstrap 95% CI excluding 0, per arm (real, sim) and pooled; else negative.
G1b: positive control - sim arm top-5 >= 90% (spike-in recovery).
G4: report-quality metrics (trace completeness, calibration, false reassurance)."""
import csv, json, os, random

BASE = "/home/sandbox/rare-disease-variant-pipeline"
RES = os.path.join(BASE, "repo/results")

def recall(rows, col, k):
    ok = [r for r in rows if r["causal_in_panel"] == "True"]
    hit = sum(1 for r in ok if r[col] and int(r[col]) <= k)
    return hit / len(ok) if ok else 0.0, hit, len(ok)

def boot_diff(rows, k, n_boot=10000, seed=3):
    rng = random.Random(seed)
    ok = [r for r in rows if r["causal_in_panel"] == "True"]
    n = len(ok)
    diffs = []
    for _ in range(n_boot):
        samp = [ok[rng.randrange(n)] for _ in range(n)]
        f = sum(1 for r in samp if r["full_rank"] and int(r["full_rank"]) <= k) / n
        b = sum(1 for r in samp if r["baseline_rank"] and int(r["baseline_rank"]) <= k) / n
        diffs.append(f - b)
    diffs.sort()
    return sum(diffs) / n_boot, diffs[int(0.025 * n_boot)], diffs[int(0.975 * n_boot)]

if __name__ == "__main__":
    rows = list(csv.DictReader(open(os.path.join(RES, "case_level_results.tsv")),
              delimiter="\t"))
    out = {"n_cases": len(rows)}
    for arm in ("real", "sim", "all"):
        sub = rows if arm == "all" else [r for r in rows if r["arm"] == arm]
        f5, fh, fn = recall(sub, "full_rank", 5)
        b5, bh, bn = recall(sub, "baseline_rank", 5)
        f1, _, _ = recall(sub, "full_rank", 1)
        b1, _, _ = recall(sub, "baseline_rank", 1)
        mean_d, lo, hi = boot_diff(sub, 5)
        out[arm] = {"n": len(sub),
                    "full_top5": round(f5, 4), "base_top5": round(b5, 4),
                    "full_top1": round(f1, 4), "base_top1": round(b1, 4),
                    "top5_diff": round(mean_d, 4), "ci95": [round(lo, 4), round(hi, 4)],
                    "g3_pass": (mean_d >= 0.05 and lo > 0)}
    # G1b positive control: sim-arm full top5 >= 0.90
    out["g1b_pass"] = out["sim"]["full_top5"] >= 0.90
    # G4: calibration + false reassurance + trace completeness
    top1_p = [r for r in rows if r["top1_class"] in ("P", "LP")]
    out["g4"] = {
        "top1_plp_rate": round(len(top1_p) / len(rows), 4),
        "false_reassurance": round(sum(1 for r in rows
            if r["causal_full_class"] in ("B", "LB")) / len(rows), 4),
        "causal_class_dist": {c: sum(1 for r in rows if r["causal_full_class"] == c)
                              for c in ("P", "LP", "VUS", "LB", "B")},
    }
    json.dump(out, open(os.path.join(RES, "g3_g4_metrics.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))
