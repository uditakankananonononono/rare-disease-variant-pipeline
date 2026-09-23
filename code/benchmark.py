#!/usr/bin/env python3
"""G1/G2 variant-level gate evaluation from results_variant_level.tsv.
G1a: tuning-split balanced accuracy >= 0.75 (P/LP-vs-rest sensitivity and
B/LB-vs-rest specificity; VUS calls count as errors on both sides).
G2: held-out MCC(full) - MCC(baseline) > 0 with bootstrap 95% CI excluding 0,
else the negative is reported."""
import csv, json, os, random, math

BASE = "/home/sandbox/rare-disease-variant-pipeline"
RES = os.path.join(BASE, "results_variant_level.tsv")

def load():
    return list(csv.DictReader(open(RES), delimiter="\t"))

def metrics(rows, arm_col):
    tp = fp = tn = fn = 0
    for r in rows:
        truth = int(r["label"])
        call = 1 if r[arm_col] in ("P", "LP") else 0
        if truth == 1 and call == 1: tp += 1
        elif truth == 1 and call == 0: fn += 1
        elif truth == 0 and call == 1: fp += 1
        else: tn += 1
    sens = tp / (tp + fn) if tp + fn else 0
    spec = tn / (tn + fp) if tn + fp else 0
    bal = (sens + spec) / 2
    denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = ((tp * tn - fp * fn) / denom) if denom else 0
    vus_rate = sum(1 for r in rows if r[arm_col] == "VUS") / len(rows)
    return {"tp": tp, "fp": fp, "tn": tn, "fn": fn, "sens": round(sens, 4),
            "spec": round(spec, 4), "balanced_acc": round(bal, 4),
            "mcc": round(mcc, 4), "vus_rate": round(vus_rate, 4)}

def bootstrap_mcc_diff(rows, n_boot=10000, seed=1):
    rng = random.Random(seed)
    n = len(rows)
    diffs = []
    for _ in range(n_boot):
        samp = [rows[rng.randrange(n)] for _ in range(n)]
        mf = metrics(samp, "full_label")["mcc"]
        mb = metrics(samp, "baseline_label")["mcc"]
        diffs.append(mf - mb)
    diffs.sort()
    lo, hi = diffs[int(0.025 * n_boot)], diffs[int(0.975 * n_boot)]
    return sum(diffs) / n_boot, lo, hi

if __name__ == "__main__":
    rows = load()
    tune = [r for r in rows if r["split"] == "tune"]
    held = [r for r in rows if r["split"] == "heldout"]
    out = {"n_total": len(rows), "n_tune": len(tune), "n_heldout": len(held)}
    for name, subset in (("tune", tune), ("heldout", held)):
        out[name] = {"baseline": metrics(subset, "baseline_label"),
                     "full": metrics(subset, "full_label")}
    mean_d, lo, hi = bootstrap_mcc_diff(held)
    out["g2_bootstrap"] = {"mcc_diff_mean": round(mean_d, 4),
                           "ci95": [round(lo, 4), round(hi, 4)],
                           "pass": lo > 0}
    out["g1a_pass"] = out["tune"]["full"]["balanced_acc"] >= 0.75
    os.makedirs(os.path.join(BASE, "repo/results"), exist_ok=True)
    json.dump(out, open(os.path.join(BASE, "repo/results/g1_g2_metrics.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))
