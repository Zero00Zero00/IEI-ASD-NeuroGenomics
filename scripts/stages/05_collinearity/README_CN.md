# MAR05 Collinearity Audit v1.0R1

Scientific purpose:
Audit overlap/collinearity among the seven raw26-derived higher-order units
separately in the PGC and SPARK mapped gene universes, and decide whether
MAR03 conditional coefficients may appear in main-text supportive reporting.

Primary metrics:
- pairwise intersection / union / Jaccard / phi
- VIF for each primary unit
- centered/scaled design-matrix condition number and singular values
- design correlation (phi)
- native coefficient-correlation availability
- secondary standalone overlap matrix

Frozen thresholds:
- VIF: GREEN <5; CAUTION 5–10; HIGH >10
- condition number: GREEN <30; CAUTION 30–100; HIGH >100
- pairwise Jaccard: GREEN <0.5; CAUTION 0.5–0.7; HIGH >0.7
- hard conditional-model rule: if any primary VIF >10 OR condition number >100,
  conditional results are SUPPLEMENT_ONLY; otherwise MAIN_TEXT_SUPPORTIVE_ELIGIBLE.
- marginal models remain PRIMARY regardless.

No association scores or beta values are used in the collinearity computations.
The runner directly calls the already validated MAR02/MAR03 Python executable;
no `conda run` overhead.

Recommended:
1. `bash run_mar05.sh all`
2. Review `MAR05_preanalysis_lock.json`
3. `bash run_mar05.sh release RELEASE`
4. `bash run_mar05.sh all-post`
5. Upload `MAR05_REVIEW_<UTC>.tar.gz` + `.sha256`
