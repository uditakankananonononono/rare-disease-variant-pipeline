# DATA MANIFEST — byte-lock (milestone 2), 2026-09-23 (Asia/Calcutta)

All raw files live outside git (sandbox data dir); this manifest + SHA-256s are the byte-lock.

## Downloaded artifacts (SHA-256)
| file | source URL | sha256 |
|---|---|---|
| variant_summary.txt.gz (443,063,763 B) | https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz | f25970aa49619fc73cdf8c878a3a20faaa69a0abd56cc8e14954f1d143e7e44e |
| genes_to_phenotype.txt | http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt | 507a17bff9c49e6329fbd88b1f91734fa7e9fdea5c06304044c0892eb6ab248c |
| phenotype.hpoa | http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa | e89aa39c8f97bf5a52c5d160f9250a603681d630ade8ec2679d84ac5aece1f72 |
| clinvar_grch38_small.tsv (derived: GRCh38, SNV/ins/del/dup/indel <=50bp, 8 germline sig classes; 4,172,171 rows) | derived from variant_summary.txt.gz (filter in code/ notes below) | ded13ae7787005b251cda304976b0244aa048d06057975b1bceab88d98ad4e88 |

## Git-pinned artifacts
| artifact | pin |
|---|---|
| phenopacket-store (monarch-initiative) | commit 4aed56ed8b1cdcf100336172a13cf144db7a5a60, cloned 2026-09-23; 10,714 phenopackets / 730 gene cohorts |
| phenopacket_case_index.tsv (derived, in this repo) | 11,086 causal-variant records; 9,738 evaluable (hg38, <=50bp, ACMG P/LP, >=3 observed HPO terms, gene mapped) |

## Remote-slice artifacts (queried at runtime; per-query log kept)
| source | endpoint / path |
|---|---|
| GIAB HG002 NIST v4.2.1 GRCh38 benchmark VCF (+ .tbi) | https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz |
| gnomAD v4 GraphQL (allele freq, gene constraint) | https://gnomad.broadinstitute.org/api |
| Ensembl REST VEP (consequence, SIFT, PolyPhen) | https://rest.ensembl.org/vep/human/region |

## ClinVar subset class composition (GRCh38 small variants)
Pathogenic 176,434; Pathogenic/Likely pathogenic 39,775; Likely pathogenic 118,078;
Uncertain significance 2,318,458; Conflicting 164,785; Likely benign 1,084,815;
Benign/Likely benign 66,182; Benign 203,644.
Review status: expert panel 21,437; practice guideline 22; multiple submitters no
conflict 664,719; single submitter 3,223,966; no assertion 97,600; conflicting 164,427.

## Known caveats
- ClinVar answer key is label, not truth; review-status >= 2 stars required for eval splits (G5).
- phenopacket ACMG labels are curator assertions; used as case-truth for G3 real arm.
- GIAB HG002 background: one real individual; population-background realism caveat stated in paper.
