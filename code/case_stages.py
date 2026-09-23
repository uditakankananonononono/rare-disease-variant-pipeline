#!/usr/bin/env python3
"""Stage driver for the case-level pipeline (cache-resume-safe).
Usage: case_stages.py coords|hg002|annotate|finish"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

BASE = "/home/sandbox/rare-disease-variant-pipeline"
DATA = os.path.join(BASE, "data")

stage = sys.argv[1]
if stage == "coords":
    from build_cases import gene_coords
    genes = set()
    cdir = os.path.join(DATA, "cases")
    for fn in sorted(os.listdir(cdir)):
        c = json.load(open(os.path.join(cdir, fn)))
        genes.update(c["panel"])
    cache_fp = os.path.join(DATA, "cache/gene_coords.json")
    have = set(json.load(open(cache_fp))) if os.path.exists(cache_fp) else set()
    print(f"coords: {len(have & genes)}/{len(genes)} cached")
    gene_coords(sorted(genes))
    have = set(json.load(open(cache_fp)))
    print(f"coords now: {len(have & genes)}/{len(genes)}")
elif stage == "hg002":
    from build_cases import gene_coords
    from run_cases import prefetch_hg002, assemble_cases
    genes = set()
    cdir = os.path.join(DATA, "cases")
    cases = [json.load(open(os.path.join(cdir, fn))) for fn in sorted(os.listdir(cdir))]
    for c in cases:
        genes.update(c["panel"])
    coords = gene_coords(sorted(genes))
    t0 = time.time()
    prefetch_hg002(cases, coords)
    print(f"hg002 pass done in {time.time()-t0:.0f}s; cached genes: "
          f"{len(os.listdir(os.path.join(DATA, 'cache/hg002')))}")
elif stage == "annotate":
    from run_cases import assemble_cases, annotate_cases
    cases = assemble_cases()
    print(f"{len(cases)} cases, candidates/case ~{sum(len(c['candidates']) for c in cases)//len(cases)}")
    annotate_cases(cases)
    print("annotate pass done")
elif stage == "finish":
    print("run run_cases.py main flow via finish_run.py")
