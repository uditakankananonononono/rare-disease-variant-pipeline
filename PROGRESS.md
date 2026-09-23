# PROGRESS (auto-updated by builder 24)
- 2026-09-23 17:05 IST: gnomAD AF cache COMPLETE 6002/6000. Gene-constraint cache 1842 (annotate_split resume-safe, ~200 genes/call). Next: finish constraint -> indexes (~1 min) -> results_variant_level.tsv -> benchmark.py G1/G2 -> milestone 3 report.
  (2 parallel workers, Ensembl throttling managed with capped backoff). Next: gnomAD +
  constraint from cache, G1/G2 metrics (code/benchmark.py), then case build + G3/G4.
- Seal-criterion note: methodological contribution = RDCB two-arm case-report benchmark
  (full arm vs InterVar-automated-tier baseline arm) + report-quality metric suite
  (trace completeness, calibration, false-reassurance) quantified on 6,000 ClinVar
  variants + 400 cases. User grant on connected accounts/internet tools noted.
