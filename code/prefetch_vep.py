#!/usr/bin/env python3
"""Parallel VEP cache prefetch. Usage: prefetch_vep.py WORKER_IDX N_WORKERS
Shards 50-variant batches across workers; cache files make merging automatic."""
import csv, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import features
from features import vep_batch

BASE = "/home/sandbox/rare-disease-variant-pipeline"
w_idx, n_w = int(sys.argv[1]), int(sys.argv[2])
rows = list(csv.DictReader(open(os.path.join(BASE, "repo/data/variant_eval_split.tsv")),
            delimiter="\t"))
variants = [(r["chrom"], r["pos"], r["ref"], r["alt"]) for r in rows]
todo = []
for c, p, r, a in variants:
    fp = os.path.join(features.CACHE, "vep", f"{c}-{p}-{r}-{a}.json")
    if not os.path.exists(fp):
        todo.append((c, p, r, a))
batches = [todo[i:i + 50] for i in range(0, len(todo), 50)]
mine = [b for i, b in enumerate(batches) if i % n_w == w_idx]
logf = os.path.join(BASE, f"prefetch_vep_{w_idx}.log")
with open(logf, "a") as fh:
    fh.write(f"{time.strftime('%H:%M:%S')} start: {len(mine)} batches of {len(batches)} todo\n")
for i, b in enumerate(mine):
    try:
        vep_batch(b)
    except Exception as e:
        with open(logf, "a") as fh:
            fh.write(f"{time.strftime('%H:%M:%S')} batch {i} error: {e}\n")
        time.sleep(5)
    if i % 10 == 0:
        with open(logf, "a") as fh:
            fh.write(f"{time.strftime('%H:%M:%S')} {i}/{len(mine)}\n")
with open(logf, "a") as fh:
    fh.write(f"{time.strftime('%H:%M:%S')} WORKER DONE\n")
