# MAR02 raw26 Tier-1 redundancy compression — v1.0R1

## 基于真实 Linux 路径/接口

本执行包来自 2026-09-09 `MAR02_INTERFACE_PROBE` 的真实接口核查，而不是 SOP 中的示意接口。

真实上游：
- `/home/h3021/chapter3/12_molecular_autism_revision/01_anchor_ingest/WP01_raw26_primary_v1/`
- `WP01_raw26_tier1.tsv`
- `WP01_raw26_tier.tsv.gz`
- `WP01_raw26_mdv3_all.tsv.gz`
- `WP01_raw26_crossstage_MDV2.tsv`
- `WP01_raw26_gene_set.tsv`
- `WP01_raw26_vs_core25.tsv.gz`

真实 legacy algorithm authority：
- `/home/h3021/chapter3/config/mdv5_domain_gate_v1_1.yaml`
- legacy MDV5 scripts 02–05（本包内为 interface probe 时的冻结快照）
- historical output `/home/h3021/chapter3/05_domains/MDV5_v1_1/`

正式 MAR02 输出：
- `/home/h3021/chapter3/12_molecular_autism_revision/02_raw26_domains/MAR02_raw26_compression_v1/`

`02_raw26_domains/MA_R02_neuro_iei_primary_v1/` 是旧 CUR-R02 临床证据项目，禁止覆盖。

## 关键适配

legacy MDV5 代码内部以 `CoreSeed_flag` 记录 rare anchor。MAR02 不修改算法，而是在 `work/` 里生成一个 **raw26 adapter universe**：仅在 legacy 代码内部把 26 个 raw26 genes 映射到这个旧字段；所有正式输出均重命名为 `raw26_anchor_flag/raw26_driver_n`，避免把 CoreSeed25 误当成 primary。

所有 Leiden/threshold/稳定性参数都从真实 `/config/mdv5_domain_gate_v1_1.yaml` 继承；preflight 强制核对：J=0.70、edge=0.20、sensitivity=.15/.25/.30、100 runs、stable >=3 & >=.70、small=2 & >=.70。seed/base、resolution、objective 等不猜测，全部继承真实冻结配置。

## 最简执行

上传解压后：

```bash
cd /home/h3021/chapter3/12_molecular_autism_revision/MAR02_RAW26_COMPRESSION_v1.0R1

tmux new -s mar02
bash run_mar02.sh all
```

第一次 `all` 会完成 preflight + preanalysis freeze，然后故意停在人工 Gate A。

检查：

`/home/h3021/chapter3/12_molecular_autism_revision/02_raw26_domains/MAR02_raw26_compression_v1/MAR02_preanalysis_lock.json`

确认只读取 WP01/raw26、gene universe、pathway membership、legacy algorithm config，没有 PGC/SPARK common-side results 后：

```bash
bash run_mar02.sh release RELEASE
bash run_mar02.sh all
```

完成后将自动生成 review package，并显示路径。

## PASS 后最重要的分支

看：

`MAR02_raw26_vs_legacy_mdv5_summary.tsv`

- `ALL_DOWNSTREAM_MEMBERSHIP exact_same=1`：MAR03 不重做 SNP harmonization/MAGMA gene mapping，只从冻结 gene-level results 重跑 gene-set 层。
- `exact_same=0`：raw26-derived units 为 primary，legacy MDV5 降为 historical sensitivity；MAR03–MAR08 按新 units 运行。
