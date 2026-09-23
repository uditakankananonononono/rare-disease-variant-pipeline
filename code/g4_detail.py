#!/usr/bin/env python3
"""G4 detail: trace completeness, calibration bins, false-reassurance CI, top-5 churn.
Codes may appear as PS1/PM5 or strength-suffixed (PVS1_Moderate)."""
import csv, json, os, re, random

BASE = "/home/sandbox/rare-disease-variant-pipeline"
RES = os.path.join(BASE, "repo/results")
CODE_RE = re.compile(r'((?:PVS1|PS|PM|PP|BA1|BS|BP)\w*)\s*:')

def main():
    rows = list(csv.DictReader(open(os.path.join(RES, "case_level_results.tsv")), delimiter="\t"))
    checked = incomplete = no_code = 0
    for r in rows:
        txt = open(os.path.join(RES, "case_reports", r["case_id"] + ".md")).read()
        line = re.search(r"^\| 1 \|.*$", txt, re.M).group(0)
        trace = line.split("|")[-2].strip()
        checked += 1
        if trace == "no criteria fired":
            no_code += 1
        elif not CODE_RE.findall(trace):
            incomplete += 1
    bins = {}
    for r in rows:
        txt = open(os.path.join(RES, "case_reports", r["case_id"] + ".md")).read()
        conf = re.search(r"^\| 1 \|.*$", txt, re.M).group(0).split("|")[-3].strip()
        hit = bool(r["full_rank"] and int(r["full_rank"]) == 1)
        b = bins.setdefault(conf, [0, 0]); b[0] += hit; b[1] += 1
    rng = random.Random(7)
    fr = [1 if r["causal_full_class"] in ("B", "LB") else 0 for r in rows]
    n = len(fr)
    boots = sorted(sum(rng.choice(fr) for _ in range(n)) / n for _ in range(10000))
    lost = sum(1 for r in rows if r["baseline_rank"] and int(r["baseline_rank"]) <= 5
               and not (r["full_rank"] and int(r["full_rank"]) <= 5))
    gained = sum(1 for r in rows if r["full_rank"] and int(r["full_rank"]) <= 5
                 and not (r["baseline_rank"] and int(r["baseline_rank"]) <= 5))
    out = {
        "trace": {"checked": checked, "complete": checked - incomplete,
                  "no_fired_codes": no_code, "completeness": (checked - incomplete) / checked},
        "calibration_top1": {k: {"hits": v[0], "n": v[1], "acc": round(v[0] / v[1], 4)}
                             for k, v in sorted(bins.items())},
        "false_reassurance": {"rate": sum(fr) / n, "ci95_boot": [boots[250], boots[9750]],
                              "rule_of_three_upper": round(3 / n, 4)},
        "top5_churn": {"lost_by_full": lost, "gained_by_full": gained},
    }
    json.dump(out, open(os.path.join(RES, "g4_detail.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))

if __name__ == "__main__":
    main()
