# PROGRESS (auto-updated by builder 24)
- 2026-09-23 15:58 IST: VEP cache 4652/6000 via foreground chunks (sandbox suspends background jobs between turns; chunk protocol in wake prompt). Ensembl ~15s/request floor; ~1300 variants left, then annotate_split (gnomAD/constraint) -> benchmark.py G1/G2 -> milestone 3 report.
  (2 parallel workers, Ensembl throttling managed with capped backoff). Next: gnomAD +
  constraint from cache, G1/G2 metrics (code/benchmark.py), then case build + G3/G4.
- Seal-criterion note: methodological contribution = RDCB two-arm case-report benchmark
  (full arm vs InterVar-automated-tier baseline arm) + report-quality metric suite
  (trace completeness, calibration, false-reassurance) quantified on 6,000 ClinVar
  variants + 400 cases. User grant on connected accounts/internet tools noted.
