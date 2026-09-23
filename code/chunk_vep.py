#!/usr/bin/env python3
"""Foreground VEP prefetch with time budget. Usage: chunk_vep.py MAX_SECONDS"""
import csv, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import features
from features import vep_batch

BASE = "/home/sandbox/rare-disease-variant-pipeline"
w_idx = int(sys.argv[1]) if len(sys.argv) > 2 else 0
n_w = int(sys.argv[2]) if len(sys.argv) > 2 else 1
budget = float(sys.argv[3]) if len(sys.argv) > 3 else 100
t0 = time.time()
rows = list(csv.DictReader(open(os.path.join(BASE, "repo/data/variant_eval_split.tsv")),
            delimiter="\t"))
variants = [(r["chrom"], r["pos"], r["ref"], r["alt"]) for r in rows]
todo = []
for c, p, r, a in variants:
    fp = os.path.join(features.CACHE, "vep", f"{c}-{p}-{r}-{a}.json")
    if not os.path.exists(fp):
        todo.append((c, p, r, a))
print(f"todo: {len(todo)}")
done = 0
batches = [todo[i:i + 50] for i in range(0, len(todo), 50)]
for bi, b in enumerate(batches):
    if bi % n_w != w_idx:
        continue
    if time.time() - t0 > budget:
        break
    try:
        vep_batch(b)
        done += len(b)
    except Exception as e:
        print(f"batch error: {e}")
        time.sleep(3)
print(f"processed {done} in {time.time()-t0:.0f}s; cache now "
      f"{len(os.listdir(os.path.join(features.CACHE, 'vep')))}/6000")
