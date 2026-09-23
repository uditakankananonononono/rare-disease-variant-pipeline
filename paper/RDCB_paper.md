# RDCB: A Case-Level Interpretable Reporting Benchmark for Rare-Disease Variant Interpretation, Evaluated on Public Undiagnosed-Cohort Data

## Abstract

Rare-disease diagnosis from genomic data is usually evaluated as variant classification
or variant ranking. What a clinician or a family actually receives, however, is neither:
it is a case-level report - a ranked, argued, evidence-linked statement about which
variant in this patient explains the disease and why. No public benchmark evaluates
that artifact. We built RDCB (Rare-Disease Case-Bench): an open pipeline and benchmark
whose evaluated artifact is the interpretable case report. Every candidate in every
report carries a per-evidence ACMG-2015 trace linking each fired criterion to a
retrievable public datum, plus an inheritance hypothesis and a calibrated confidence
statement. The benchmark runs entirely on public data: 200 real diagnosed cases from
the GA4GH phenopacket-store and 200 simulated cases spiking graded ClinVar pathogenic
variants into the GIAB HG002 v4.2.1 high-confidence variant background, evaluated over
disease-gene panels derived from each case's HPO terms (2,648 unique genes, 21,440
unique candidate variants). All gates were locked in writing before any outcome data
were seen, and both negative results are reported unaltered. Headline gate G3 failed:
adding phenotype/inheritance integration to a faithful InterVar-style ACMG baseline
did not improve top-5 causal-variant recall (pooled difference -1.82 percentage
points, 95% CI [-3.66, -0.26]). The pre-registered secondary metric showed the
opposite at rank 1: +13.4 percentage points top-1 recall on real cases (0.8656 vs
0.7312). The phenotype term sharpens the single best answer and slightly blurs the
shortlist - a quantified design tradeoff for case-ranking systems. Report quality,
the benchmark's other axis, held: 100% evidence-trace completeness across 400 case
reports, top-1 confidence calibration 1.00/1.00/0.845 across high/moderate-high/
uncertain bins, and a false-reassurance rate of 0/400 (< 0.75% upper bound).
Pipeline, data manifests with SHA-256 checksums, all metrics, and all 400 case
reports are public.

## 1. Introduction: the evaluated artifact is wrong

Genomic rare-disease diagnosis still ends without an answer for most families; cohort
solves cluster around 25-50% depending on phenotype and platform. The field's
evaluation culture tracks two artifacts: per-variant classification (ClinVar-style
label prediction) and per-case ranking (is the causal variant near the top of a
list). Neither is what is handed to a clinician. A diagnostic report is an argument:
it must say which variant, in which gene, with which evidence, under which
inheritance model, at what confidence - and every step of that argument must be
checkable against the underlying datum. Automation that cannot show its evidence is
not auditable; a ranking metric cannot tell whether the report a family receives is
trustworthy.

We asked a simple question: can an open, fully reproducible pipeline, using only
public data and public APIs, produce case-level diagnostic reports whose evidence is
completely traceable - and does integrating phenotype and inheritance information
actually help find the causal variant, measured honestly against a stated baseline
under gates locked before scoring?

The answer is mixed, and both halves matter. The headline comparison failed: at
top-5, phenotype integration did not help. The secondary comparison succeeded: at
top-1, it helped a great deal. And the report-quality artifact - the thing no prior
benchmark measures - passed its gates. We report all of it, including two negative
results, exactly as the pre-locked gates require.

## 2. Prior art: crowded at classification and ranking, empty at reports

The rare-disease interpretation space is crowded (full analysis in PRIOR_ART.md;
verdict CROWDED with the case-report angle surviving). MOLGENIS VIP (NAR Genomics &
Bioinformatics 2025, doi.org/10.1093/nargab/lqaf087) automates per-variant
decision-tree classification and emits a variant-table HTML report, but is not
benchmarked at case-report level. AutoGVP (Kim et al., Bioinformatics 2024,
doi.org/10.1093/bioinformatics/btae114) automates ClinVar+InterVar ACMG calls per
variant. Exomiser (Robinson et al., 2020, PMC7230372) and LIRICAL
(doi.org/10.1101/2020.01.25.19014803) perform phenotype-driven prioritization;
LIRICAL's likelihood ratios are disease-level, not ACMG evidence-traced per variant.
VarPPUD (doi.org/10.1002/humu.24459) benchmarks prioritization; PhEval
(s12859-025-06105-4) benchmarks rankers. nf-core/raredisease ships WGS/WES workflows.
None evaluates the interpretable case report - ranked candidates with per-evidence
traceability, inheritance hypothesis, and calibrated confidence - as the benchmarked
artifact on public undiagnosed-cohort data. That artifact is RDCB's contribution.

## 3. Data: multiple public sources, byte-locked

All inputs are public and were byte-locked with SHA-256 checksums before any scoring
(DATA_MANIFEST.md); the manifest was re-verified after a full sandbox restore:

- ClinVar variant_summary.txt.gz (graded germline classifications; used ONLY as
  answer key for variant-level gates and same-codon evidence for OTHER variants -
  never as a feature for the evaluated variant itself).
- gnomAD v4 API: allele frequencies and gene constraint (LOEUF, mis_z).
- Ensembl VEP REST: consequence, SIFT, PolyPhen, transcript identifiers.
- HPO: genes_to_phenotype.txt and phenotype.hpoa (case phenotype to gene mapping).
- GIAB HG002 v4.2.1 benchmark VCF (GRCh38, chromosomes 1-22): the real individual
  variant background for simulated cases.
- GA4GH phenopacket-store (monarch-initiative, commit 4aed56e): 11,086 curated
  published case records; 9,738 evaluable with HPO terms and a causal variant.

Two cohorts were built. Real arm: 200 diagnosed phenopacket cases (186 with the
causal gene inside the derived panel; the 14 panel-miss cases are retained as honest
unsolvable cases and reported separately). Sim arm: 200 simulated cases, each
spiking one graded (2-star+) ClinVar P/LP variant into the HG002 background within
the case's HPO-derived panel (196 causal-in-panel). Panels derive from each case's
HPO terms via genes_to_phenotype; across 400 cases the panels cover 2,648 unique
genes and 21,440 unique candidate variants (mean 61 candidates/case), each annotated
with VEP consequence, gnomAD frequency, and gene constraint (100% coverage;
per-variant cached, failures never fabricated).

Variant-level evaluation uses an independent stratified split: 6,000 ClinVar
variants (2-star+ review, Mendelian genes, balanced P/LP vs B/LB, conflicts
excluded), 4,742 tuning / 1,258 held-out, split by variant and documented.

Known data limitations, stated in advance: GRCh38 SNVs and small indels (<=50 bp)
only; no SVs/CNVs/repeats; GIAB v4.2.1 GRCh38 covers autosomes 1-22 only, so
simulated-case backgrounds carry no chrX variants (148 panel genes affected,
recorded per gene in the cache manifest); gnomAD constraint is population-
ascertainment-limited for recessive and late-onset genes (Section 5.2).

## 4. Methods

### 4.1 Two arms, gates locked first

The baseline arm is our faithful implementation of InterVar's published automated
ACMG-2015 rule tier (automated criteria only, no manual adjustment) - a stated
single-tool baseline chosen in writing before scoring. The full arm keeps the
baseline's variant classification and adds case-level integration: phenotype-gene
relevance (HPO term-set overlap between the case and each candidate gene),
inheritance/zygosity fit, and gnomAD constraint carried as REPORTED EVIDENCE in the
trace - not as downgrade criteria (Amendment 1, Section 5.2). The evaluated artifact
is the case report: candidates ranked with per-evidence traces.

Gates were locked in GATES_LOCKED.md (committed before any results commit, git
history public): G1 positive control (variant-level balanced accuracy >= 0.75 on the
tuning split; causal variant top-5 in >= 90% of spike-in controls), G2 (held-out
variant-level full-vs-baseline MCC, bootstrap 95% CI excludes 0), G3 (headline:
full-arm top-5 causal recall must beat baseline by >= 5 absolute points with
bootstrap CI excluding 0, on >= 200 real + >= 200 sim cases; top-1 recall is the
pre-registered secondary metric), G4 (report quality: 100% top-1 evidence-trace
completeness; calibration and false-reassurance reported descriptively), G5 (no
leakage: no feature derives from the evaluated variant's own ClinVar record), G6
(honesty: splits, exclusions, and verified/thin/missing counts published; negatives
preserved, never re-fished).

### 4.2 Case construction and ranking

Real cases use the phenopacket's HPO terms and causal variant directly. Simulated
cases take a real phenopacket's HPO profile and spike its graded ClinVar causal
variant into HG002's PASS calls within the derived panel, so the background is a
real individual's variation, not simulated genotypes. Baseline-arm ranking orders by
ACMG class, then allele frequency. Full-arm ranking adds the phenotype-overlap score
to the class score. Each case report lists every candidate with class, confidence,
and an evidence trace in which every fired ACMG code is linked to its datum
(e.g., "PM2: gnomAD max AF 0.000007 < 5e-5"; "PVS1_Moderate: LoF in LoF-mechanism
gene AIRE but LOEUF upper 1.055 >= 0.35"), plus a provenance block naming every
feature source.

## 5. Results

### 5.1 G1 (positive control): pass

Baseline-arm variant-level balanced accuracy on the tuning split: 0.8124 (sensitivity
0.6249, specificity 1.0, MCC 0.6742, VUS rate 0.5803) - above the 0.75 gate.
End-to-end spike-in control (G1b): full-arm top-5 causal recovery 0.9796 on 200
simulated cases, above the 0.90 gate (baseline arm 1.0).

### 5.2 G2 (held-out variant-level): negative, preserved

The original full arm (v1) modulated PVS1 by gene constraint: PVS1 fired at full
strength only when the gene's LOEUF upper bound was < 0.35, else at Moderate
strength. On the held-out split this destroyed sensitivity without buying anything:
baseline MCC 0.678 (sensitivity 0.6302, specificity 1.0) vs full-arm v1 MCC 0.285
(sensitivity 0.1508); MCC difference -0.393, bootstrap 95% CI [-0.431, -0.355]. Gate
G2 failed; the negative is reported as-is (AMENDMENT_v1_negative.md, written before
any case-level run).

Diagnosis: of 1,197 true P/LP variants the full arm downgraded on the tuning split,
1,194 were the LOEUF demotion. Recessive and late-onset disease genes are
LoF-TOLERANT in population data precisely because LoF carriers are healthy - yet LoF
is the disease mechanism. gnomAD constraint, designed to find genes intolerant to
variation, is the wrong tool for demoting LoF evidence in Mendelian disease genes.
This misuse pattern is tempting (constraint looks like a mechanistic prior) and we
have not seen the failure quantified elsewhere; we report it as a design warning.

Per Amendment 1, the case-level full arm was then redefined (baseline classification
+ phenotype/inheritance integration + constraint as reported evidence only), and G2
stands failed against v1.

### 5.3 G3 (headline, case-level): negative, preserved

| Arm | n (evaluable) | baseline top-5 | full top-5 | diff (pp) | 95% CI | G3 |
|---|---|---|---|---|---|---|
| real | 200 (186) | 0.9892 | 0.9731 | -1.61 | [-4.30, +1.08] | fail |
| sim | 200 (196) | 1.0000 | 0.9796 | -2.04 | [-4.08, -0.51] | fail |
| pooled | 400 (382) | 0.9948 | 0.9764 | -1.82 | [-3.66, -0.26] | fail |

Pre-registered secondary metric, top-1 recall: full 0.8656 vs baseline 0.7312 on real
cases (+13.4 pp); 0.9337 vs 0.8724 on sim (+6.1 pp); pooled +9.7 pp.

Churn analysis: across 382 evaluable cases, the full arm's ranking pushed the causal
variant out of the top-5 in 12 cases and pulled it in in 2. The phenotype-overlap
term promotes the true causal variant to rank 1 - its effect on the top of the list
is real and large - but it also lifts phenotype-matched non-causal VUSs into ranks
2-5, and causal variants in genes with weak or missing HPO-to-gene mappings receive
no phenotype score and sink. Net: better single best answer, slightly worse
shortlist. For a report meant to be read, rank 1 carries special weight; for a
shortlist meant to be worked through, top-5 is the metric. The pre-locked headline
gate chose top-5, and by that gate the integration failed. Both numbers stand.

### 5.4 G4 (report quality): gates met

Across all 400 case reports:

- Evidence-trace completeness: 400/400 = 100%. Every fired ACMG code on every top-1
  report links to its retrievable datum (verified programmatically, including
  strength-suffixed codes such as PVS1_Moderate with the LOEUF value shown).
- Confidence calibration (stated top-1 confidence vs empirical top-1 accuracy):
  high 1/1 = 1.000; moderate-high 77/77 = 1.000; uncertain 272/322 = 0.845.
  Stated confidence is conservative: when the report says high or moderate-high, it
  was always right in this cohort.
- False-reassurance rate (B/LB call on the true causal variant): 0/400 = 0.0
  (upper bound < 0.75% by the rule of three).
- Causal-variant class distribution under the full arm: P 1, LP 77, VUS 322, LB 0,
  B 0. The pipeline re-derives classification from evidence rather than importing
  ClinVar labels (G5), so most known-causal variants surface as VUS with their
  partial evidence shown - itself an honest picture of what public evidence
  supports variant-by-variant.

## 6. The tool

RDCB ships as an open pipeline (code/) plus the benchmark artifacts: 400 case
reports (results/case_reports/), case-level and variant-level result tables, all
gate metrics as JSON, the locked gates, the amendment, and byte-locked data
manifests. Reports are Markdown, human-readable by construction, and each names its
feature provenance. The pipeline is cache-resumable end to end (per-variant VEP and
gnomAD caches; per-gene GIAB slices; constraint cache) and was restored from the
manifest once after a full environment wipe with all checksums re-verified.

## 7. Discussion

Three findings survive contact with the locked gates.

First, the interpretable case report is a benchmarkable artifact, and it can be
produced at 100% evidence-trace completeness from public data alone. Prior
benchmarks score rankings; RDCB scores the argument. Calibration results say the
argument's confidence language means something: stated high and moderate-high
confidences were never wrong in this cohort.

Second, phenotype integration is a top-1 intervention, not a top-5 one. The measured
tradeoff (+13.4 pp top-1, -1.6 pp top-5 on real cases) is a concrete design result
for case-ranking systems: if the product is a single best answer with a trace,
integrate phenotype; if the product is a shortlist for manual review, rank on
variant evidence and show phenotype as context. RDCB's full arm does the former and
its reports show both.

Third, gene constraint is dangerous as a downgrade criterion for Mendelian LoF
variants (G2): the population-genetics meaning of LOEUF is nearly opposite to the
mechanistic question being asked. Kept as reported evidence (Amendment 1), it is
informative without being destructive.

Limitations. SNV/small-indel scope only; no chrX in simulated backgrounds (GIAB
v4.2.1 GRCh38 is autosomal); real-arm diagnoses come from published case reports,
which skew toward solved, exome-amenable disease; the pipeline's VUS-heavy class
distribution reflects strict evidence re-derivation without clinical correlate data;
Ensembl/gnomAD API drift means exact numbers are tied to the byte-locked manifests
and cached evidence, all of which ship with the benchmark.

## 8. Reproducibility

Repository: code/, data/ (manifests + SHA-256 checksums), results/, paper/.
GATES_LOCKED.md predates every results commit (git history). RESUME.md documents
full environment restore; one restore was performed and all checksums re-verified.
All external data are public; no wet-lab work; no proprietary features.

## 9. References

1. MOLGENIS VIP. NAR Genom Bioinform 2025. doi.org/10.1093/nargab/lqaf087
2. Kim et al. AutoGVP. Bioinformatics 2024. doi.org/10.1093/bioinformatics/btae114
3. Robinson et al. LIRICAL. doi.org/10.1101/2020.01.25.19014803
4. Exomiser. Robinson et al. 2020. PMC7230372
5. VarPPUD benchmark. doi.org/10.1002/humu.24459
6. PhEval. BMC Bioinformatics 2025. s12859-025-06105-4
7. InterVar automated ACMG interpretation. Wang et al. PMC7257575
8. ClinVar variant_summary. ncbi.nlm.nih.gov/clinvar
9. gnomAD v4. gnomad.broadinstitute.org
10. Ensembl VEP REST. rest.ensembl.org
11. Human Phenotype Ontology. hpo.jax.org
12. GIAB HG002 v4.2.1 benchmark. ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab
13. GA4GH phenopacket-store. github.com/monarch-initiative/phenopacket-store
14. nf-core/raredisease. nf-co.re/raredisease
