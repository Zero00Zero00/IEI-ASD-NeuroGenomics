# MAR08_FAILURE_LOCALIZATION_AND_INFLUENCE_v1.0R1.1_AUTHORITY_FIX

本包只修复 MAR08 R1 的 upstream canonical source-data 路径假设错误；冻结的
`MAR08_REVISED_ANALYSIS_CONTRACT_v1.0` **没有改变**。

## R1 HOLD 根因

R1 假设以下 canonical TSV 已存在于 active stage 的 `source_data/`：
- MAR02 5 个 canonical source TSV
- MAR03 primary higher-order source TSV
- MAR04 writeup summary source TSV

真实 Linux 中这些旧阶段的 checksum ledger 存在，但上述 canonical source TSV 是在对应
FREEZE_RELEASE 中冻结，而非 materialized 到 active stage 目录。因此 R1 在 scientific analysis
开始之前正确 HOLD。

R1.1 将这 7 个文件作为 byte-identical frozen transport snapshots 随包携带，并：
1. 对每个文件做 exact SHA256；
2. 对来源 MAR02/MAR03/MAR04 FREEZE_RELEASE 外层 SHA256 与 frozen contract 做 provenance 对照；
3. 仍检查 active upstream stage 的 checksum ledger 存在；
4. 不重跑、不改写任何 upstream stage；
5. 使用新 output namespace `MAR08_failure_localization_influence_v1p1`。

## 执行

```bash
cd /home/h3021/chapter3/12_molecular_autism_revision
sha256sum -c MAR08_FAILURE_LOCALIZATION_AND_INFLUENCE_v1.0R1.1_AUTHORITY_FIX.tar.gz.sha256
tar -xzf MAR08_FAILURE_LOCALIZATION_AND_INFLUENCE_v1.0R1.1_AUTHORITY_FIX.tar.gz
cd MAR08_FAILURE_LOCALIZATION_AND_INFLUENCE_v1.0R1.1_AUTHORITY_FIX
sha256sum -c PACKAGE_CHECKSUMS.sha256

bash run_mar08_r1p1.sh all
```

第一次应停在 Gate A：

```text
MAR08_PREFLIGHT=PASS
SCENARIOS=L1:11 L2:36
MAR07_MDV3_SIDECARS=1000/1000
MAR08_FREEZE=PASS
HOLD_FOR_MANUAL_GATE
```

然后检查：
- `MAR08_preanalysis_lock.json`
- `MAR08_FROZEN_SCENARIO_MANIFEST.tsv`
- `MAR08_FROZEN_SCENARIO_COVERAGE.tsv`
- `MAR08_preflight_report.tsv`

人工确认后：

```bash
bash run_mar08_r1p1.sh release RELEASE
bash run_mar08_r1p1.sh all-post
```

最终上传：
`90_handoff/MAR08_REVIEW_<UTC>.tar.gz` 与 `.sha256`。
