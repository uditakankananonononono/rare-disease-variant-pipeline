#!/usr/bin/env python3
"""G3/G4 case-level benchmark: annotate candidates (shared disk cache),
classify both arms, rank, score causal-variant recall, generate per-case
interpretable reports with full evidence traces."""
import csv, json, os, sys, time, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from features import vep_batch, gnomad_batch, constraint_batch, build_codon_index, same_codon, vid
from acmg import classify, build_lof_genes, SCORE
import pysam

BASE = "/home/sandbox/rare-disease-variant-pipeline"
DATA = os.path.join(BASE, "data")
RES = os.path.join(BASE, "repo/results")
os.makedirs(RES, exist_ok=True)
os.makedirs(os.path.join(RES, "case_reports"), exist_ok=True)

def log(m):
    print(time.strftime("%H:%M:%S"), m, flush=True)
    with open(os.path.join(BASE, "cases.log"), "a") as fh:
        fh.write(time.strftime("%H:%M:%S ") + m + "\n")

# ---------- per-gene HG002 background cache ----------
GIAB_URL = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
import os as _os
_LOCAL = _os.path.join(DATA, "cache/HG002.vcf.gz")
GIAB_SRC = _LOCAL if _os.path.exists(_LOCAL) else GIAB_URL

def prefetch_hg002(cases, coords):
    d = os.path.join(DATA, "cache/hg002")
    os.makedirs(d, exist_ok=True)
    genes = sorted({g for c in cases for g in c["panel"]})
    import pysam
    tb = None
    done = 0
    for g in genes:
        fp = os.path.join(d, g + ".json")
        if os.path.exists(fp):
            done += 1
            continue
        c = coords.get(g)
        if not c or not c.get("start"):
            json.dump([], open(fp, "w")); continue
        if tb is None:
            tb = pysam.TabixFile(GIAB_SRC)
        if ("chr" + c["chrom"]) not in tb.contigs:
            # definitive: GIAB v4.2.1 GRCh38 benchmark covers autosomes 1-22 only (no chrX)
            json.dump([], open(fp, "w")); continue
        out = None
        for attempt in range(3):
            try:
                if attempt:
                    tb = pysam.TabixFile(GIAB_SRC)
                cur = []
                for row in tb.fetch("chr" + c["chrom"], int(c["start"]), int(c["stop"])):
                    f = row.split("\t")
                    if f[6] != "PASS":
                        continue
                    ref, alts = f[3], f[4].split(",")
                    for alt in alts:
                        if alt in (".", "*") or max(len(ref), len(alt)) > 50:
                            continue
                        cur.append({"gene": g, "chrom": c["chrom"], "pos": int(f[1]),
                                    "ref": ref, "alt": alt, "origin": "HG002"})
                out = cur
                break
            except Exception:
                import time as _t; _t.sleep(1)
        if out is None:
            continue
        json.dump(out, open(fp, "w"))
        done += 1
        if done % 50 == 0:
            log(f"hg002 prefetch {done}/{len(genes)}")
    log(f"hg002 prefetch complete {done}/{len(genes)}")

def load_term2genes():
    import csv as _csv
    t2g = collections.defaultdict(set)
    with open(os.path.join(DATA, "genes_to_phenotype.txt")) as fh:
        for row in _csv.DictReader(fh, delimiter="\t"):
            t2g[row["hpo_id"]].add(row["gene_symbol"])
    return t2g

def assemble_cases():
    cdir = os.path.join(DATA, "cases")
    cases = []
    for fn in sorted(os.listdir(cdir)):
        if fn.endswith(".json"):
            cases.append(json.load(open(os.path.join(cdir, fn))))
    # attach background from per-gene cache (replacing build-time fetch)
    bgcache = os.path.join(DATA, "cache/hg002")
    for c in cases:
        bg = []
        excl = {(v["chrom"], v["pos"]) for v in c["causal"]}
        for g in c["panel"]:
            fp = os.path.join(bgcache, g + ".json")
            if os.path.exists(fp):
                for v in json.load(open(fp)):
                    if (v["chrom"], v["pos"]) not in excl:
                        bg.append(v)
        import random
        random.Random(42).shuffle(bg)
        c["candidates"] = bg[:60] + [dict(v, origin="causal") for v in c["causal"]]
    return cases

def annotate_cases(cases):
    allv = {}
    for c in cases:
        for v in c["candidates"]:
            allv[vid(str(v["chrom"]), v["pos"], v["ref"], v["alt"])] = \
                (str(v["chrom"]), str(v["pos"]), v["ref"], v["alt"])
    variants = list(allv.values())
    log(f"unique candidate variants across cases: {len(variants)}")
    t0 = time.time()
    vep = vep_batch(variants); log(f"vep {time.time()-t0:.0f}s")
    gn = gnomad_batch(variants); log(f"gnomad {time.time()-t0:.0f}s")
    genes = sorted({v.get("gene", "") for v in vep.values() if v.get("gene")})
    con = constraint_batch(genes); log(f"constraint {len(genes)} genes")
    return vep, gn, con

def rank_case(c, feats, t2g):
    """Returns list of ranked candidate dicts for each arm."""
    case_hpos = set(c["hpo"])
    rows = []
    for v in c["candidates"]:
        k = vid(str(v["chrom"]), v["pos"], v["ref"], v["alt"])
        f = feats.get(k, {})
        gname = f.get("gene", "") or v.get("gene", "")
        n_overlap = sum(1 for h in case_hpos if gname in t2g.get(h, ()))
        pheno = n_overlap / max(1, len(case_hpos))
        rows.append({"vid": k, "gene": gname, "origin": v.get("origin", ""),
                     "af": f.get("af") if f.get("af") is not None else 0.0,
                     "feat": f, "pheno": pheno})
    base = sorted(rows, key=lambda r: (-SCORE[r["feat"]["baseline_label"]],
                                       r["af"], r["vid"]))
    full = sorted(rows, key=lambda r: (-(SCORE[r["feat"]["full_label"]] + r["pheno"]),
                                       r["af"], r["vid"]))
    return base, full

def report_case(c, full_rank, t2g):
    """Markdown case report: ranked candidates with per-evidence trace."""
    lines = [f"# Case report: {c['case_id']}", "",
             f"- Arm: {c['arm']} | Disease: {c['disease']} | Causal gene: {c['causal_gene']}",
             f"- HPO terms ({len(c['hpo'])}): {', '.join(c['hpo'][:12])}" +
             (" ..." if len(c['hpo']) > 12 else ""),
             f"- Panel: {len(c['panel'])} genes | Candidates evaluated: {len(c['candidates'])}",
             "", "## Ranked candidates (full arm)", "",
             "| Rank | Gene | Variant | Class | Confidence | Evidence trace |",
             "|---|---|---|---|---|---|"]
    conf_map = {"P": "high", "LP": "moderate-high", "VUS": "uncertain",
                "LB": "low (likely not causal)", "B": "very low"}
    for i, r in enumerate(full_rank[:10], 1):
        f = r["feat"]
        ev = "; ".join(f"{k}: {v}" for k, v in (f.get("full_codes") or {}).items()) or "no criteria fired"
        lines.append(f"| {i} | {r['gene']} | {r['vid']} | {f['full_label']} | "
                     f"{conf_map[f['full_label']]} | {ev} |")
    lines += ["", "## Provenance",
              "Features: Ensembl VEP REST (consequence/SIFT/PolyPhen), gnomAD v4 "
              "(AF, constraint), ClinVar same-codon evidence (other records only). "
              "Classification: two-arm ACMG-2015 automated rules (see GATES_LOCKED.md)."]
    path = os.path.join(RES, "case_reports", c["case_id"] + ".md")
    open(path, "w").write("\n".join(lines))
    return path

def main():
    from build_cases import gene_coords
    cases = assemble_cases() if os.path.exists(os.path.join(DATA, "cases")) else []
    if not cases:
        log("no cases found - run build_cases.py first")
        sys.exit(1)
    log(f"{len(cases)} cases loaded")
    coords = gene_coords(sorted({g for c in cases for g in c["panel"]}))
    prefetch_hg002(cases, coords)
    cases = assemble_cases()  # re-assemble now cache is complete
    vep, gn, con = annotate_cases(cases)
    ci = build_codon_index()
    lg = build_lof_genes()
    t2g = load_term2genes()

    feats = {}
    for c in cases:
        for v in c["candidates"]:
            k = vid(str(v["chrom"]), v["pos"], v["ref"], v["alt"])
            if k in feats:
                continue
            vv = vep.get(k, {}); gg = gn.get(k, {"af": None})
            cg = con.get(vv.get("gene", ""), {})
            ps1, pm5 = same_codon(ci, vv.get("gene", ""), vv.get("protein_start"),
                                  vv.get("amino_acids", ""), k)
            feat = {**vv, **gg, **cg, "ps1": ps1, "pm5": pm5}
            bl, bc = classify(feat, lg, "baseline")
            fl, fc = classify(feat, lg, "full")
            feat["baseline_label"], feat["baseline_codes"] = bl, bc
            feat["full_label"], feat["full_codes"] = fl, fc
            feats[k] = feat
    log("classification complete")

    results = []
    for c in cases:
        base_rank, full_rank = rank_case(c, feats, t2g)
        causal_vids = {vid(str(v["chrom"]), v["pos"], v["ref"], v["alt"]) for v in c["causal"]}
        def best_rank(ranked):
            for i, r in enumerate(ranked, 1):
                if r["vid"] in causal_vids:
                    return i
            return None
        br, fr = best_rank(base_rank), best_rank(full_rank)
        report_case(c, full_rank, t2g)
        top1 = full_rank[0]
        results.append({"case_id": c["case_id"], "arm": c["arm"],
                        "causal_gene": c["causal_gene"],
                        "causal_in_panel": c["causal_in_panel"],
                        "n_candidates": len(c["candidates"]),
                        "baseline_rank": br, "full_rank": fr,
                        "top1_class": top1["feat"]["full_label"],
                        "causal_full_class": next((feats[v]["full_label"] for v in causal_vids if v in feats), "")})
    with open(os.path.join(RES, "case_level_results.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(results[0].keys()))
        w.writeheader(); w.writerows(results)
    log(f"DONE -> {RES}/case_level_results.tsv")

if __name__ == "__main__":
    main()
