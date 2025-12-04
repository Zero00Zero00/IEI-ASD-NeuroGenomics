# 05_gene2func_axis_validation

Tissue expression and pathway validation for CoreSeed, Set A and axis gene sets.

This module validates gene sets and axes using tissue expression and independent pathway analyses.

Steps:
- Build gene sets for CoreSeed, Set A and axis-specific clusters.
- Construct expression matrices from brain and immune tissues (e.g. GTEx) and perform differential expression and tissue enrichment analyses.
- Run additional pathway enrichment analyses using these gene sets.
- Map enriched pathways back onto the three mechanistic axes to check consistency with the bias-aware enrichment module.

The outputs support the biological interpretation of the axes and are used in figures showing tissue and pathway-level validation.
