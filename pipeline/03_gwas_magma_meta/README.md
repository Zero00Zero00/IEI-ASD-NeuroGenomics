# 03_gwas_magma_meta

QC and MAGMA gene-based meta-analysis of ASD GWAS summary statistics.

This module runs QC and MAGMA-based meta-analysis of ASD GWAS summary statistics.

Steps:
- Harmonise and filter ASD GWAS summary statistics, aligning to HapMap3 SNPs and a common reference panel.
- Perform MAGMA gene-based tests across different TSS-centred windows (e.g. 0 kb, 10 kb, 50 kb).
- Implement a correlation-aware meta-analysis across cohorts using the strategy described in the manuscript.
- Export gene-level association tables for use in functional annotation and Set A construction.

These outputs form the gene-based GWAS backbone of the study and are referenced by the SNP2GENE, Set A and S-LDSC modules.
