# 04_fuma_snp2gene

FUMA-style SNP2GENE mapping integrating positional, eQTL and 3D chromatin evidence.

This module implements a FUMA-style SNP2GENE mapping pipeline.

Steps:
- Standardise CoreSeed and candidate gene lists and prepare long-format eQTL tables for multiple tissues.
- Map GWAS SNPs to genes through positional proximity, eQTL and 3D chromatin (Hi-C) annotations.
- Integrate these evidence sources into gene-level scores and tables.
- Generate intermediate and final files that summarise SNP-to-gene links and support mechanistic interpretation.

The resulting SNP2GENE annotations are used to populate the multi-modal evidence matrix and to refine Set A membership.
