#!/usr/bin/env python3
"""Build benchmark cases (G3). Real arm: phenopacket-store published cases.
Sim arm: ClinVar 2-star+ P/LP variants (disjoint from the variant-level split,
G5) spiked into GIAB HG002 background. Panel genes from case HPO terms only
(no leakage of the causal gene beyond what HPO implies). Background: HG002
NIST v4.2.1 PASS variants in panel-gene regions via remote tabix."""
import csv, json, os, random, sys, collections
import pysam
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from features import http_post

BASE = "/home/sandbox/rare-disease-variant-pipeline"
DATA = os.path.join(BASE, "data")
GIAB = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
COORD_CACHE = os.path.join(DATA, "cache/gene_coords.json")
MAX_PANEL = 40
MAX_BG = 150

INH_DOM, INH_REC = "HP:0000006", "HP:0000007"
SKIP_HPO = {INH_DOM, INH_REC, "HP:0000005"}  # inheritance terms are not phenotypes

def load_maps():
    term2genes = collections.defaultdict(set)
    with open(os.path.join(DATA, "genes_to_phenotype.txt")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            term2genes[row["hpo_id"]].add(row["gene_symbol"])
    dis2hpo, dis2inh = collections.defaultdict(list), collections.defaultdict(set)
    with open(os.path.join(DATA, "phenotype.hpoa")) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            dis, hpo = f[0], f[3]
            if hpo in (INH_DOM, INH_REC):
                dis2inh[dis].add(hpo)
            else:
                dis2hpo[dis].append(hpo)
    return term2genes, dis2hpo, dis2inh

def gene_coords(genes):
    cache = json.load(open(COORD_CACHE)) if os.path.exists(COORD_CACHE) else {}
    todo = [g for g in genes if g and g not in cache]
    for i in range(0, len(todo), 20):
        chunk = todo[i:i + 20]
        q = "query { " + " ".join(
            f'g{j}: gene(gene_symbol: "{g}", reference_genome: GRCh38) '
            "{ symbol chrom start stop }" for j, g in enumerate(chunk)) + " }"
        res = http_post("https://gnomad.broadinstitute.org/api", {"query": q})
        for j, g in enumerate(chunk):
            d = res.get("data", {}).get(f"g{j}")
            cache[g] = {"chrom": (d.get("chrom") or "").replace("chr", ""),
                        "start": d.get("start"), "stop": d.get("stop")} if d else None
        import time; time.sleep(0.5)
    json.dump(cache, open(COORD_CACHE, "w"))
    return cache

_tb = None
def tabix():
    global _tb
    if _tb is None:
        _tb = pysam.TabixFile(GIAB)
    return _tb

def background_variants(coords, panel, exclude_positions):
    bg = []
    for g in panel:
        c = coords.get(g)
        if not c or not c.get("start"):
            continue
        try:
            rows = tabix().fetch("chr" + c["chrom"], int(c["start"]), int(c["stop"]))
        except Exception:
            continue
        for row in rows:
            f = row.split("\t")
            if f[6] != "PASS":
                continue
            pos = int(f[1])
            ref, alts = f[3], f[4].split(",")
            for alt in alts:
                if alt in (".", "*") or max(len(ref), len(alt)) > 50:
                    continue
                if (c["chrom"], pos) in exclude_positions:
                    continue
                bg.append({"gene": g, "chrom": c["chrom"], "pos": pos,
                           "ref": ref, "alt": alt, "origin": "HG002"})
    random.Random(42).shuffle(bg)
    return bg[:MAX_BG]

def make_panel(case_hpos, term2genes):
    scored = collections.Counter()
    for h in case_hpos:
        for g in term2genes.get(h, ()):
            scored[g] += 1
    return [g for g, _ in scored.most_common(MAX_PANEL)]

def build_real(n_target=200):
    term2genes, dis2hpo, dis2inh = load_maps()
    idx = list(csv.DictReader(open(os.path.join(DATA,
        "../repo/data/phenopacket_case_index.tsv")), delimiter="\t"))
    idx = [r for r in idx if r["evaluable"] == "1"]
    by_case = collections.defaultdict(list)
    for r in idx:
        by_case[r["phenopacket_id"]].append(r)
    cands = [v for v in by_case.values() if int(v[0]["n_hpo"]) >= 5]
    random.Random(7).shuffle(cands)
    per_gene = collections.Counter()
    cases = []
    for variants in cands:
        gene0 = variants[0]["cohort"]  # cohort dir name is the gene symbol
        if per_gene[gene0] >= 2:
            continue
        hpos = [h for h in variants[0]["hpo_ids"].split(";") if h not in SKIP_HPO]
        panel = make_panel(hpos, term2genes)
        if len(panel) < 3:
            continue
        per_gene[gene0] += 1
        cases.append({"case_id": variants[0]["phenopacket_id"], "arm": "real",
                      "hpo": hpos, "panel": panel,
                      "causal": [{"gene": v["gene"], "chrom": v["chrom"],
                                  "pos": int(v["pos"]), "ref": v["ref"],
                                  "alt": v["alt"], "zygosity": v["allelic_state"],
                                  "hgvs_c": v["hgvs_c"]} for v in variants],
                      "causal_gene": gene0, "disease": variants[0]["disease_id"]})
        if len(cases) >= n_target:
            break
    return cases

def build_sim(n_target=200):
    term2genes, dis2hpo, dis2inh = load_maps()
    split_vids = set()
    with open(os.path.join(BASE, "repo/data/variant_eval_split.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            split_vids.add("-".join([r["chrom"], r["pos"], r["ref"], r["alt"]]))
    GOOD = {"reviewed by expert panel", "practice guideline",
            "criteria provided, multiple submitters, no conflicts"}
    pool = []
    with open(os.path.join(DATA, "clinvar_grch38_small.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["ClinicalSignificance"] not in ("Pathogenic", "Likely pathogenic"):
                continue
            if r["ReviewStatus"] not in GOOD:
                continue
            omim = next((p.split(":")[1] for p in r["PhenotypeIDS"].split(",")
                         if p.startswith("OMIM:") and "OMIM:" + p.split(":")[1] in dis2hpo), None)
            if not omim:
                continue
            vid = "-".join([r["Chromosome"], r["PositionVCF"], r["ReferenceAlleleVCF"],
                            r["AlternateAlleleVCF"]])
            if vid in split_vids or "na" in vid:
                continue
            pool.append({"gene": r["GeneSymbol"], "chrom": r["Chromosome"],
                         "pos": int(r["PositionVCF"]), "ref": r["ReferenceAlleleVCF"],
                         "alt": r["AlternateAlleleVCF"], "omim": "OMIM:" + omim,
                         "name": r["Name"]})
    random.Random(11).shuffle(pool)
    seen_gene = set()
    cases = []
    for v in pool:
        if v["gene"] in seen_gene:
            continue
        terms = [t for t in dis2hpo[v["omim"]]]
        if len(terms) < 3:
            continue
        panel = make_panel(terms, term2genes)
        if len(panel) < 3:
            continue
        inh = dis2inh.get(v["omim"], set())
        zyg = "homozygous" if (INH_REC in inh and INH_DOM not in inh) else "heterozygous"
        seen_gene.add(v["gene"])
        cases.append({"case_id": f"sim_{v['gene']}_{v['pos']}", "arm": "sim",
                      "hpo": terms, "panel": panel,
                      "causal": [{"gene": v["gene"], "chrom": v["chrom"],
                                  "pos": v["pos"], "ref": v["ref"], "alt": v["alt"],
                                  "zygosity": zyg, "hgvs_c": v["name"]}],
                      "causal_gene": v["gene"], "disease": v["omim"]})
        if len(cases) >= n_target:
            break
    return cases

if __name__ == "__main__":
    real = build_real(200)
    sim = build_sim(200)
    outdir = os.path.join(DATA, "cases")
    os.makedirs(outdir, exist_ok=True)
    coords_needed = {g for c in real + sim for g in c["panel"]}
    print("fetching coords for", len(coords_needed), "panel genes")
    coords = gene_coords(sorted(coords_needed))
    n_with_bg = 0
    with open(os.path.join(DATA, "cases_index.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["case_id", "arm", "causal_gene", "disease", "n_hpo",
                    "panel_size", "causal_in_panel", "n_background", "n_candidates"])
        for c in real + sim:
            excl = {(v["chrom"], v["pos"]) for v in c["causal"]}
            bg = background_variants(coords, c["panel"], excl)
            cands = bg + [dict(v, origin="causal") for v in c["causal"]]
            c["candidates"] = cands
            c["causal_in_panel"] = c["causal_gene"] in c["panel"]
            if bg:
                n_with_bg += 1
            json.dump(c, open(os.path.join(outdir, c["case_id"] + ".json"), "w"))
            w.writerow([c["case_id"], c["arm"], c["causal_gene"], c["disease"],
                        len(c["hpo"]), len(c["panel"]), c["causal_in_panel"],
                        len(bg), len(cands)])
    print("real:", len(real), "sim:", len(sim), "with background:", n_with_bg)
    real_in = sum(1 for c in real if c["causal_in_panel"])
    sim_in = sum(1 for c in sim if c["causal_in_panel"])
    print(f"causal gene in panel: real {real_in}/{len(real)}, sim {sim_in}/{len(sim)}")
