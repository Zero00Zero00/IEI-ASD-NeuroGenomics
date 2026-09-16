#!/usr/bin/env python3
from common import *
ensure_dirs()
if not pass_file(WORK/'MAR02_PREFLIGHT_PASS.txt'): raise SystemExit('HOLD: run preflight first')
config=WORK/'MAR02_adapter_config.yaml'
files=[WP01/'WP01_raw26_gene_set.tsv',WP01/'WP01_raw26_tier1.tsv',WP01/'WP01_raw26_tier.tsv.gz',WP01/'WP01_raw26_mdv3_all.tsv.gz',WP01/'WP01_raw26_crossstage_MDV2.tsv',WP01/'WP01_raw26_vs_core25.tsv.gz',UNIVERSE,MEMBERSHIP,LEGACY_CFG,config,WORK/'MAR02_raw26_universe_adapter.tsv',WORK/'MAR02_raw26_tier1_adapter.tsv']
files += sorted(VENDOR.glob('mdv5_*.py'))
rules={'redundancy_complete_linkage_jaccard':0.70,'representative_rule':'q_emp -> q_Firth -> pathway_size -> source_priority -> pathway_id','primary_edge_threshold':0.20,'sensitivity_edge_thresholds':[0.15,0.25,0.30],'leiden_runs':100,'stable_domain':'>=3 representatives and median coclustering >=0.70','small_module':'exactly 2 representatives and pair coclustering >=0.70','Expanded':'union of representative pathway protein-coding genes','Compact':'genes in >=2 representative pathways; sensitivity only','primary_common_family':'higher-order Expanded units','secondary_family':'standalone representatives','magma_gene_window_kb':10}
obj={'stage':'MAR02','version':VERSION,'timestamp_utc':utc(),'association_results_read':'NO','gwas_inputs_used':'NO','common_side_files_read':'NO','output_root':str(OUT),'rules':rules,'input_hashes':{str(p):sha256_file(p) for p in files if p.is_file()},'forbidden_common_side_inputs':['/06_magma_pgc/','/07_magma_spark/','PGC2019_gene_results_w10.tsv','SPARK_EUR_gene_results_w10.tsv']}
atomic_json(OUT/'MAR02_preanalysis_lock.json',obj)
atomic_text(OUT/'MAR02_FREEZE_READY.txt',f'status=PASS\nstage=MAR02_FREEZE_READY\npreanalysis_lock_sha256={sha256_file(OUT/"MAR02_preanalysis_lock.json")}\nmanual_gate_required=YES')
print('MAR02_FREEZE=PASS\nMANUAL_GATE_REQUIRED=YES\nNEXT=bash run_mar02.sh release RELEASE')
