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
