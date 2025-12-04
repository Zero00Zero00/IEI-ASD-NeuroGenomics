# 02_pathway_bias_and_CPS11

Bias-aware pathway enrichment and CPS-11 pathway prioritisation.

This module performs bias-aware pathway enrichment and defines the mechanistic axes.

Steps:
- Construct the protein-coding gene universe and annotate genes with length, GC content and other structural features.
- Run multiple pathway enrichment frameworks (standard ORA, logistic regression with splines, Mantel–Haenszel, goseq) to correct for gene-length and related biases.
- Cluster enriched GO/Reactome terms into pathway communities and assign them to candidate axes.
- Compute CPS-11 scores to prioritise robust pathways and summarise them at the axis level.

The outputs provide the pathway-level evidence and axis definitions used in later modules (Set A, axis validation and S-LDSC).
