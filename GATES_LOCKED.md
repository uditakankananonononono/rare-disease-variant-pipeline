# RDCB — Rare-Disease Case-Bench: LOCKED GATES (written 2026-09-23, before any outcome data)

## Contribution (named before gates lock)
An open benchmark + pipeline where the evaluated artifact is the **case-level
interpretable diagnostic report** — not a variant classifier. Each case report ranks
candidate variants with per-evidence ACMG-2015 code traceability (every code linked to a
retrievable datum), an inheritance hypothesis, and a calibrated confidence statement.
The benchmark runs on **public undiagnosed-cohort data**: (a) published case reports
curated as GA4GH phenopackets (monarch-initiative/phenopacket-store), each with HPO terms
+ causal variant, and (b) simulated cases spiking graded ClinVar pathogenic variants into
a real individual's variant background (GIAB HG002 v4.2.1 high-confidence calls in
HPO-derived panel genes). Diagnostic-accuracy metrics (top-k recall of the causal
variant) AND report-quality metrics (evidence-trace completeness, confidence calibration,
false-reassurance rate) are reported together. Report quality has never been the
benchmarked artifact in this field (see PRIOR_ART.md).

## Prior-art verdict: CROWDED (angle verified surviving)
Closest works: MOLGENIS VIP (per-variant decision-tree classification + variant-table
HTML report; no case-report benchmark), AutoGVP (per-variant ClinVar+InterVar ACMG
automation; no case-level artifact), Exomiser/LIRICAL (phenotype-driven ranking;
LIRICAL's LR output is disease-level, not an ACMG evidence-traced case report),
VarPPUD (ranking), PhEval (benchmarks rankings, not reports). None ships a benchmark
whose evaluated artifact is the interpretable case report on public undiagnosed-cohort
data. Full citations in PRIOR_ART.md.

## Scope lock
- Variant scope: SNVs + small indels (<=50bp), GRCh38, nuclear genome, Mendelian disease
  genes with an OMIM/HPO disease mapping. No SVs/CNVs/repeats (stated limitation).
- Feature sources (all public): ClinVar variant_summary (graded labels = answer key
  ONLY, never a feature for the variant itself), gnomAD v4 API (allele frequencies,
  gene constraint LOEUF), Ensembl VEP REST (consequence, SIFT, PolyPhen), HPO
  (genes_to_phenotype, phenotype.hpoa), GIAB HG002 v4.2.1 benchmark VCF (background).
- AlphaMissense/EVE: used only if stream-filterable within RAM; features, never truth.

## Baseline (stated in writing before scoring)
Single-tool baseline arm = our faithful implementation of InterVar's published automated
ACMG rule tier (Wang et al., PMC7257575 logic: automated criteria only, no manual
adjustment). Full arm adds gene constraint (PP2/BP1/PVS1-decision support), same-codon
ClinVar evidence (PS1/PM5), phenotype-gene relevance (HPO term-set overlap with IC
weighting), and inheritance-model fit at case level.

## Gates
- G1 (positive control, tuning split only): (a) variant-level balanced accuracy >= 0.75
  vs ClinVar labels on the tuning split; (b) end-to-end: causal variant ranked top-5 in
  >= 90% of positive-control spike-in cases. Failing G1 = pipeline bug, fix and re-run;
  G1 numbers reported either way.
- G2 (held-out variant-level): on the held-out split (stratified random 20%, review
  status >= 2 stars, P/LP vs B/LB, conflicts excluded), full-arm MCC must exceed
  baseline-arm MCC with a bootstrap 95% CI on the difference excluding 0, else the
  negative is reported.
- G3 (held-out case-level, headline): on >= 200 evaluable phenopacket-store cases plus
  >= 200 simulated spike-in cases, full-arm causal-variant top-5 recall must exceed the
  baseline arm by >= 5 absolute percentage points (bootstrap 95% CI excludes 0), else
  negative. Top-1 recall reported as secondary.
- G4 (report quality): 100% of top-1 reports must carry a complete evidence trace (each
  fired ACMG code linked to its datum); confidence calibration reported (binned stated
  vs empirical accuracy); false-reassurance rate (B/LB call on a truly pathogenic causal
  variant) reported with 95% CI — no performance claim made on G4, it is descriptive.
- G5 (leakage): no feature derives from the evaluated variant's own ClinVar record;
  same-codon evidence uses other variants only; split is by variant, documented.
- G6 (honesty): all splits, exclusion counts, and verified/thin/missing evidence counts
  published in every report. Negatives preserved, never re-fished.

## Data byte-lock plan (milestone 2)
ClinVar variant_summary.txt.gz (release dated), phenotype.hpoa, genes_to_phenotype.txt,
phenopacket-store git commit SHA, GIAB HG002 VCF byte-range manifest (per-gene tabix
slices actually read), gnomAD/Ensembl API query logs. Each local artifact: SHA-256.
