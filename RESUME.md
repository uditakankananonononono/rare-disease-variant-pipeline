# RESUME / fresh-sandbox restore (builder 24)

1. Re-download raw data per data/DATA_MANIFEST.md and verify SHA-256 (ClinVar
   variant_summary.txt.gz, genes_to_phenotype.txt, phenotype.hpoa; phenopacket-store
   at commit 4aed56e) into /home/sandbox/rare-disease-variant-pipeline/data/.
2. Rebuild ClinVar subset: filter = Assembly==GRCh38, Type in
   {single nucleotide variant, Deletion, Insertion, Duplication, Indel}, len<=50bp,
   8 germline significance classes; expect 4,172,171 rows, sha256
   ded13ae7787005b251cda304976b0244aa048d06057975b1bceab88d98ad4e88.
3. pip install biopython pandas numpy matplotlib scipy pysam
4. Rebuild case index: python3 code/parse_phenopackets.py <pps/notebooks> data/phenopacket_case_index.tsv
5. Rebuild split: python3 code/build_split.py (deterministic, md5-seeded)
6. Feature annotation resumes from data/cache: python3 code/annotate_split.py
   (idempotent; caches VEP/gnomAD/constraint per-variant JSON)
7. Cases: python3 code/build_cases.py (deterministic seeds 7/11/42)
8. Benchmark + report generator: code/benchmark.py, code/report.py (pending)
Status at last commit: annotate_split running (background, cached); build_cases
smoke-tested (real arm 53/60 causal-in-panel, median panel 40 genes).
