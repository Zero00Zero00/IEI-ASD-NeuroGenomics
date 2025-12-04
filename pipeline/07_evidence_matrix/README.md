# 07_evidence_matrix

Construction of the multimodal evidence matrix for Set A and axis-level summaries.

This module builds the multi-modal evidence matrix for Set A and summarises axis-level support.

Steps:
- Merge GWAS gene-level statistics, CoreSeed flags, mechanistic axes and SNP2GENE-based positional / eQTL / 3D evidence.
- Derive composite evidence categories and scores for each Set A gene.
- Summarise evidence by axis and by mechanistic class.
- Produce heatmaps and summary plots corresponding to the Set A evidence matrix figures in the manuscript.

The resulting tables in `annotations/` can be used directly to prioritise genes or to construct axis-aware gene panels and scores.
