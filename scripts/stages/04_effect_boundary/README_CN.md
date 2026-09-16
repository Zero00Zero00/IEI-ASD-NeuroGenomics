# MAR04 Effect Boundary v1.0R1

输入：冻结的 `MAR03_FREEZE_RELEASE_v1.0R2`

固定公式：
- 95% CI = beta ± z(0.975) × SE
- one-sided 95% upper bound = beta + z(0.95) × SE
- 80% MDE = [z(1-alpha_family)+z(0.80)] × SE
- alpha_family = 0.05/m (one-sided Bonferroni)
- Primary higher-order family m=7
- Secondary standalone family m=13
- No prespecified SESOI; equivalence claims are prohibited.
- Conditional results remain SUPPORTIVE_PENDING_MAR05.

本包只使用 Python 标准库，不需要 conda/pandas/scipy。

推荐：
1. `bash run_mar04.sh all`
2. 它在 freeze 后停止。
3. 检查 `MAR04_preanalysis_lock.json`
4. `bash run_mar04.sh release RELEASE`
5. `bash run_mar04.sh all-post`

输出目录：
`/home/h3021/chapter3/12_molecular_autism_revision/04_effect_boundary/MAR04_effect_boundary_v1`

最终 `source_data/` 可直接用于 forest plot、Results 与 MAR05。
