# Prior-art check — 2026-09-23 — Verdict: CROWDED (named angle survives)

## Named tools checked (amended-bar requirement)
1. **MOLGENIS VIP** — https://molgenis.github.io/vip/ ; NAR Genom Bioinform 2025,
   doi.org/10.1093/nargab/lqaf087. End-to-end interpretation pipeline: VEP annotation,
   CAPICE prioritization, customizable decision-tree variant classification, stand-alone
   HTML report (variant table + genome browser). The report is a per-variant browser with
   per-variant classifications — NOT a case-level diagnostic narrative, and VIP is not
   benchmarked with the case report as the evaluated artifact on undiagnosed-cohort data.
2. **nf-core/raredisease** — https://nf-co.re/raredisease/3.1.2/docs/output/. WGS/WES
   calling + VEP/CADD annotation + GENMOD ranking. Output: ranked VCFs + QC. No
   interpretable case report artifact, no case-report benchmark.
3. **AutoGVP** — Kim et al., Bioinformatics 2024, doi.org/10.1093/bioinformatics/btae114.
   Dockerized per-variant ACMG classification (ClinVar + InterVar + AutoPVS1). Per-variant
   only; no case-level report; benchmark is variant-classification concordance.

## Adjacent works (closest to the surviving angle)
4. **LIRICAL** — Robinson et al., doi.org/10.1101/2020.01.25.19014803 (Am J Hum Genet
   2020). Interpretable likelihood-ratio CASE-LEVEL output — closest existing work — but
   it prioritizes DISEASES from phenotypes (+optional genotype LR), does not produce an
   ACMG evidence-traced variant interpretation report, and its benchmarks measure
   diagnostic recall, not report quality.
5. **Exomiser** — Smedley et al., Nat Protoc 2015; benchmarked on real patient WES
   (Robinson et al. 2020, PMC7230372). Ranked variant/gene lists with scores; HTML
   reports are result tables, not evidence-traced case narratives.
6. **VarPPUD** — PLoS Comput Biol 2025, journal.pcbi.1013414. Pinpoints diagnostic
   variants from prioritized candidate sets. A ranking method, benchmarked on ranking.
7. **PhEval** — BMC Bioinformatics 2025, s12859-025-06105-4. Standard benchmark harness
   for phenotype-driven gene/variant prioritization — evaluates rankings, not reports.
8. **genomicsITER/benchmark-germline-variants-prioritizers** (Hum Mutat 2022,
   doi.org/10.1002/humu.24459) — benchmark of prioritizers; ranking metrics only.
9. Commercial (Congenica, Franklin/Genoox, Fabric GEM) ship case-level reports but are
   closed and publish no open case-report benchmark on public undiagnosed cohorts.

## Conclusion
The named three tools ship per-variant classification/ranking (+VIP's variant-table
report). No open tool ships a benchmark whose evaluated artifact is the case-level
interpretable report (evidence-traced ACMG reasoning + inheritance + calibrated
confidence) scored on public undiagnosed-cohort data (phenopacket-store published cases +
GIAB-background spike-ins). Angle survives. BUILD (no KILL).
