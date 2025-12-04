# IEI-ASD-NeuroGenomics

A reproducible framework for prioritizing neuro-immune risk in autism spectrum disorder using bias-decoupled genetics.

---

## 1. Overview

This repository contains the code and curated resources for the study:

> **“A Bias-Decoupled Genetic Scaffold Links Inborn Errors of Immunity to the Polygenic Architecture of Autism”**

In this work we use inborn errors of immunity (IEI) as a natural experiment to define a compact set of *IEI–ASD CoreSeed* genes, derive three convergent neuro-immune axes, and show that common-variant ASD risk is disproportionately concentrated in the genomic neighbourhoods of these axes. We then organise the resulting gene sets into an axis-anchored resource that can be reused for risk stratification, pathway-aware analyses and trial enrichment in ASD and related neuro-immune conditions.

The repository is organised to make the full analysis **transparent and reusable**:

- A modular **Snakemake pipeline** that mirrors the Methods, from CoreSeed construction and bias-aware pathway analysis to ASD GWAS meta-analysis, functional annotation and stratified LD-score regression.
- Curated **annotation tables** for the IEI–ASD CoreSeed, Set A / Set B, and the three mechanistic axes.
- Reviewer-friendly **notebooks** that reproduce the main figures and summary statistics from static TSV files.

> Primary GWAS summary statistics, LD reference panels and expression resources (e.g. GTEx, PsychENCODE) are **not** redistributed here for licensing reasons. We instead provide download links and version numbers so that users can obtain these data from the original sources.

---

## 2. Repository layout

The repository is intended to stay lightweight. Only code, configuration and small annotation files are tracked; large raw data and derived outputs are ignored.

```text
IEI-ASD-NeuroGenomics/
├── README.md
├── LICENSE                      # to be added (e.g. MIT/Apache 2.0)
├── environment.yml              # base conda environment
├── pipeline/                    # modular Snakemake workflows
│   ├── 01_coreseed_and_mechanism/
│   ├── 02_pathway_bias_and_CPS11/
│   ├── 03_gwas_magma_meta/
│   ├── 04_fuma_snp2gene/
│   ├── 05_gene2func_axis_validation/
│   ├── 06_setA_setB_construction/
│   ├── 07_evidence_matrix/
│   └── 08_sldsc_partitioning/
├── annotations/                 # static resources and gene sets
│   ├── coreseed/
│   ├── pathways/
│   ├── axes/
│   ├── gwas_and_sets/
│   ├── fuma/
│   └── sldsc/
└── notebooks/                   # reviewer-oriented Jupyter notebooks
    ├── 00_QuickStart.ipynb
    ├── 01_CoreSeed_and_bias.ipynb
    ├── 02_Pathways_and_axes.ipynb
    ├── 03_SetA_multimodal_evidence.ipynb
    ├── 04_SLDSc_partitioning.ipynb
    └── 05_Translational_axes_and_targets.ipynb
