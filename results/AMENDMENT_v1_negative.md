# AMENDMENT 1 (2026-09-23, written AFTER G1/G2, BEFORE any G3 case-level run)

## What happened (honest negative, preserved)
Full-arm v1 (LOEUF-modulated PVS1: PVS1 only when LOEUF upper < 0.35, else
PVS1_Moderate) was evaluated against the locked gates:
- G1a (tuning balanced accuracy >= 0.75): baseline arm PASS (0.8124); full arm v1 FAIL (0.5730).
- G2 (held-out full-arm MCC must exceed baseline arm, bootstrap 95% CI excludes 0):
  FAIL. Baseline MCC 0.678 vs full-arm v1 MCC 0.285; diff -0.393, CI [-0.431, -0.355].
Cause (tuning-split diagnosis, 1,197 downgraded true P/LP variants): 1,194 are the
LOEUF>=0.35 PVS1 demotion; 3 missing-LOEUF. Recessive/late-onset disease genes are
LoF-TOLERANT in gnomAD yet LoF is the disease mechanism, so the demotion destroys
sensitivity (0.630 -> 0.151) while specificity stays 1.0. Full-arm v1 gains: 62.

## Decision
The G2 negative is reported as-is; full-arm v1 is NOT re-fished. For the G3 case-level
benchmark (not yet run; no case-level outcomes seen), the full arm is redefined as:
  baseline-arm variant classification (the tier that passed G1a)
  + phenotype-gene relevance (case HPO overlap) + inheritance/zygosity fit
  + gnomAD constraint (LOEUF, mis_z) carried as REPORTED EVIDENCE in the case-report
    trace, not as classification-downgrade criteria.
The methodological claim shifts accordingly: the evaluated contribution is case-level
integration and report quality, not an improved variant classifier. This is the
framing the amended bar specified (interpretable case-level reporting, not another
variant classifier). G3/G4 gates unchanged.
