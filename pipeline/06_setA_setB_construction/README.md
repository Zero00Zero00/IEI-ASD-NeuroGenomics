# 06_setA_setB_construction

Rules to construct the extended ASD candidate networks (Set A / Set B).

This module defines the extended ASD candidate networks (Set A and, optionally, Set B).

Steps:
- Combine MAGMA meta-analysis results with CoreSeed membership and functional annotations.
- Apply the decision rules described in the manuscript to classify genes into Set A tiers (e.g. A++, A+, A, CoreSeed-only).
- Optionally generate a broader Set B using more permissive thresholds or evidence combinations.
- Export tidy tables summarising gene IDs, tiers, axis assignments and basic features.

These sets are the main gene panels used for downstream evidence integration, S-LDSC annotation and translational mapping.
