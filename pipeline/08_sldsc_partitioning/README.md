# 08_sldsc_partitioning

Creation of SNP annotations and stratified LD-score regression analyses.

This module performs stratified LD-score regression using axis-anchored annotations.

Steps:
- Convert gene sets (CoreSeed, Set A and matched random sets) into TSS-centred genomic windows and BED files.
- Build SNP-level annotation files compatible with the baseline-LD v2.2 model across different window sizes.
- Compute LD scores for each annotation using the 1000 Genomes reference panel.
- Run S-LDSC to estimate heritability enrichment and perform sensitivity analyses, including window-size variation and leave-one-gene-out tests.
- Generate summary tables and plots for enrichment estimates and null distributions.

These analyses quantify how much ASD SNP-heritability is concentrated near axis-anchored gene sets relative to matched genomic backgrounds.
