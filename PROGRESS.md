# PROGRESS (auto-updated by builder 24)
- 2026-09-23 17:05 IST: gnomAD AF cache COMPLETE 6002/6000. Gene-constraint cache 1842 (annotate_split resume-safe, ~200 genes/call). Next: finish constraint -> indexes (~1 min) -> results_variant_level.tsv -> benchmark.py G1/G2 -> milestone 3 report.
  (2 parallel workers, Ensembl throttling managed with capped backoff). Next: gnomAD +
  constraint from cache, G1/G2 metrics (code/benchmark.py), then case build + G3/G4.
- Seal-criterion note: methodological contribution = RDCB two-arm case-report benchmark
  (full arm vs InterVar-automated-tier baseline arm) + report-quality metric suite
  (trace completeness, calibration, false-reassurance) quantified on 6,000 ClinVar
  variants + 400 cases. User grant on connected accounts/internet tools noted.

## 2026-09-23 17:28 IST
- coords stage COMPLETE: 2648/2648 panel genes resolved via gnomAD (1 unresolved symbol, skipped downstream). Incremental-save patch added after sandbox timeout killed end-only writes.
- Next: hg002 per-gene GIAB tabix prefetch (~2.6k genes), then annotate (long pole), then finish + G3/G4.

## 2026-09-23 17:39 IST
- hg002 prefetch IN PROGRESS: 1223/2648 genes cached with real tabix fetches.
- BUG FOUND+FIXED: a dropped remote tabix connection made the old except-pass path write 1520 empty cache files in 11s (false done). Deleted all dead-pass files; prefetch now retries 3x with a fresh TabixFile and leaves failures UNCACHED. Data-integrity note: dead-pass empties never reached any case or result.

## 2026-09-23 17:52 IST
- hg002 stage COMPLETE: 2648/2648 panel genes cached. Downloaded GIAB HG002 v4.2.1 VCF+tbi locally (156MB) after remote-tabix throttling; local fetch ~1ms/gene.
- Two integrity bugs fixed: (1) dead-connection pass wrote 1520 false-empty caches - purged, retry+fresh-handle patch; (2) chrX fetches always raised (GIAB v4.2.1 GRCh38 is autosomes 1-22 only, verified dir listing) - chrX now resolves to documented definitive empty, no retry. Empty caches: 148 chrX + 123 autosomal no-benchmark-call genes + 1 unresolved symbol. Limitation recorded for paper: no chrX HG002 background.
- Next: annotate stage (VEP+gnomAD for unique case candidates) - the long pole.

## 2026-09-23 18:08 IST
- annotate stage IN PROGRESS: vep cache 14902 (6002 variant-level + case candidates), gnomad 6002 pending. 21440 unique case candidates total (~61/case x 400 cases).
- Throughput fix: VEP now 6-threaded x 100/batch (5.7x: ~1500 variants/105s), gnomAD 3-threaded x 20/batch; failed batches skip uncached for resume (no fabricated annotations).

## 2026-09-23 18:42 IST
- annotate stage: vep cache 24652 (~2.8k case candidates remaining), gnomad 6002 (case pass starts when vep drains). Ensembl transient throttle resolved by dropping to 4x50 threads; 0-skip steady ~800/105s.

## 2026-09-23 19:13 IST
- annotate: vep COMPLETE for all 21,440 case candidates (cache 27,428 total). gnomad 10,662/27,428 (case pass running, ~16.8k left). gnomad throttles hard on concurrency - settled on 1 worker + 1s pacing, ~400 variants/105s, ~2 skips/call.

## 2026-09-23 19:44 IST
- annotate: gnomad 16,622/27,428 (~10.8k left, steady ~450 variants/105s, 2-6 skips/call retried next pass). Est. 2 more wake cycles to drain, then finish stage (classify/rank/report + g3g4) and M4.

## 2026-09-23 20:15 IST
- annotate: gnomad 22,622/27,428 (~4.8k left). Next wake should drain gnomad, then run finish stage (classify/rank/report + g3g4) and report M4 with G3/G4 numbers.

## 2026-09-23 20:46 IST
- annotate: vep 27,428 + gnomad 27,428 COMPLETE for all case candidates. constraint pass running (2,792 genes cached, ~240/call). Next: annotate pass done -> finish stage (classify/rank/report + g3g4) -> M4.

## 2026-09-23 20:56 IST - MILESTONE: case-level gates computed
- annotate COMPLETE (21,440 candidates: vep 0 missing, gnomad 0 missing, constraint 2,040/2,040 genes).
- run_cases DONE: 400 cases, 400 case reports with per-candidate evidence traces.
- G3 (headline, as locked): FAIL - honest negative #2. top-5 recall full vs baseline: real 0.9731 vs 0.9892 (diff -1.61pp, CI [-4.3,+1.08]); sim 0.9796 vs 1.0 (-2.04pp, CI [-4.08,-0.51]); pooled -1.82pp (CI [-3.66,-0.26]). Pre-registered secondary top-1: full +13.4pp real (0.8656 vs 0.7312), +6.1pp sim, +9.7pp pooled.
- G1b (positive control): PASS - sim full top-5 0.9796 >= 0.90.
- G4: trace completeness 400/400 (6 strength-suffixed codes were parser artifacts, each code+datum present); calibration top-1: high 1/1, moderate-high 77/77, uncertain 272/322=0.845; false reassurance 0/400 (rule-of-three upper 0.75%); causal class dist P1/LP77/VUS322/LB0/B0.
- Diagnosis: 12 lost vs 2 gained at top-5. Phenotype score boosts the true causal at rank 1 but injects pheno-matched VUS noise at ranks 2-5; causal variants in HPO-unmapped genes get pheno=0 and sink.
