# MAR06 Matched-null Audit v1.0R2 REAL

R1 HOLD 原因不是 matched-null 数据缺失，而是自动 discovery 无法区分：
- K=50 primary candidate pool vs K=20 sensitivity candidate pool；
- 10,000 primary assignments vs K20 assignments vs 1,000 pilot assignments。

上传的 `MAR06_RAW26_ARTIFACT_PROBE_20260909T140704Z` 已经明确真实冻结 authority，因此 R2 禁止启发式猜文件，改为精确绑定：

Primary K=50:
- WP01_raw26_candidate_pools_K50.tsv.gz
- WP01_raw26_full_assignments_10k.tsv.gz
- WP01_raw26_full_sets_10k.tsv.gz
- WP01_raw26_full_balance_10k.tsv.gz
- WP01_raw26_full_balance_summary.tsv

Sensitivity K=20:
- WP01_raw26_candidate_pools_K20.tsv.gz
- WP01_raw26_K20_assignments_10k.tsv.gz
- WP01_raw26_K20_sets_10k.tsv.gz
- WP01_raw26_K20_balance_10k.tsv.gz
- WP01_raw26_K20_balance_summary.tsv

Primary K50 expected from the probe:
- 26 anchors × 50 candidates = 1,300 candidate rows
- 10,000 replicates × 26 assignments = 260,000 rows
- 10,000/10,000 unique set hashes
- 0 within-set duplicate genes
- assignment-derived gene sets match frozen set files 10,000/10,000
- all 50 candidates are used for every anchor
- max anchor-level reuse proportion ~0.0256
- inverse-Simpson effective support ~49.67–49.83
- balance gates pass for all 3 structural covariates.

K20 sensitivity:
- 26 × 20 candidates
- 10,000 unique sets
- all 20 candidates used for each anchor
- max reuse proportion ~0.0571
- balance gates pass.

R2 audits existing frozen objects only; no matched set is regenerated or retuned.
Reuse concentration is diagnostic. K20 remains a prespecified sensitivity and cannot replace K50 primary.

Run:
1. `bash run_mar06.sh all`
2. Expected: authority/preflight/freeze PASS, then HOLD_FOR_MANUAL_GATE.
3. Review `MAR06_preanalysis_lock.json`.
4. `bash run_mar06.sh release RELEASE`
5. `bash run_mar06.sh all-post`
6. Upload generated `MAR06_REVIEW_<UTC>.tar.gz` + `.sha256`.

No conda-run overhead; validated Python executable is called directly.
