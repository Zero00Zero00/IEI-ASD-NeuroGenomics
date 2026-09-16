#!/usr/bin/env python3
from common import *
import platform
for p in [OUT/'MAR02_ANALYSIS_COMPLETE.txt',OUT/'MAR02_CANONICALIZE_PASS.txt',OUT/'MAR02_LEGACY_COMPARE_PASS.txt',OUT/'MAR02_SENSITIVITY_PASS.txt']:
    if not pass_file(p): raise SystemExit(f'HOLD: prerequisite missing {p}')
checks=[]
def ck(cid,ok,obs='',exp='',detail=''): checks.append([cid,'PASS' if ok else 'HOLD',str(obs),str(exp),str(detail)])
cfg=load_generated_config(); n=int(cfg['expected']['tier1_n'])
rmap=read_tsv(OUT/'MAR02_raw26_redundancy_map.tsv'); reps=read_tsv(OUT/'MAR02_raw26_representatives.tsv'); dm=read_tsv(OUT/'MAR02_raw26_domain_membership.tsv'); units=read_tsv(OUT/'MAR02_raw26_units.tsv'); st=read_tsv(OUT/'MAR02_raw26_standalones.tsv'); ptu=read_tsv(OUT/'MAR02_raw26_pathway_to_unit.tsv'); sens=read_tsv(OUT/'MAR02_raw26_threshold_sensitivity.tsv')
ck('tier1_mapped_once',len(rmap)==n and rmap.pathway_key.nunique()==n,len(rmap),n)
ck('one_rep_per_family',bool((rmap.groupby('redundancy_family_id')['representative_flag'].sum()==1).all()),int(rmap.representative_flag.sum()),rmap.redundancy_family_id.nunique())
ck('complete_linkage_minJ',bool((pd.to_numeric(rmap.family_min_pairwise_jaccard)>=0.70-1e-12).all()),rmap.family_min_pairwise_jaccard.min(),'>=0.70')
ck('pathway_to_unit_once',len(ptu)==n and ptu.pathway_key.nunique()==n,len(ptu),n)
ck('thresholds',sorted(round(float(x),12) for x in sens.edge_threshold.unique())==[0.15,0.20,0.25,0.30],sorted(sens.edge_threshold.unique()),'[.15,.20,.25,.30]')
if len(dm):
    ck('unique_unit_gene',not dm.duplicated(['unit_id','HGNC_id']).any(),int(dm.duplicated(['unit_id','HGNC_id']).sum()),0)
    ck('compact_subset_expanded',bool((pd.to_numeric(dm.loc[pd.to_numeric(dm.Compact_flag)==1,'Expanded_flag'])==1).all()),'',True)
    ck('raw26_anchor_binary',set(pd.to_numeric(dm.raw26_anchor_flag,errors='coerce').dropna().unique()).issubset({0,1}),set(pd.to_numeric(dm.raw26_anchor_flag,errors='coerce').dropna().unique()),'{0,1}')
else: ck('empty_domain_membership',len(units)==0,len(units),0)
if len(units):
    stable=units[units.unit_type=='stable_domain']; small=units[units.unit_type=='small_module']
    ck('stable_rule',bool(((pd.to_numeric(stable.representative_count)>=3)&(pd.to_numeric(stable.median_coclustering)>=0.70-1e-12)).all()),len(stable),'all pass')
    ck('small_rule',bool(((pd.to_numeric(small.representative_count)==2)&(pd.to_numeric(small.median_coclustering)>=0.70-1e-12)).all()),len(small),'all pass')
# All representatives must be either assigned to primary unit or present in standalone/unstable table.
repkeys=set(reps.pathway_key.astype(str)); primary=set(ptu.loc[ptu.unit_id.fillna('').astype(str).str.len()>0,'representative_pathway_key'].astype(str)); standalone=set(st.pathway_key.astype(str)); ck('representatives_accounted',repkeys==primary|standalone,len(repkeys-(primary|standalone)),0,';'.join(sorted(repkeys-(primary|standalone))))
# Common-side data firewall: canonical outputs/config may not carry PGC/SPARK identifiers.
forbidden=['PGC2019_gene_results_w10','SPARK_EUR_gene_results_w10','/06_magma_pgc/','/07_magma_spark/']; hits=[]
# Only configured input values are relevant; the legacy config may legitimately contain forbidden-marker strings as guardrails.
for k,v in cfg.get('upstream',{}).items():
    vv=str(v).replace('\\','/')
    for f in forbidden:
        if f.lower() in vv.lower(): hits.append(f'upstream.{k}:{f}')
for p in sorted(OUT.glob('MAR02_raw26_*')):
    if p.is_file() and p.stat().st_size<50*1024*1024:
        txt=(gzip.open(p,'rt',encoding='utf-8',errors='replace').read(2000000) if p.suffix=='.gz' else p.read_text(encoding='utf-8',errors='replace'))
        for f in forbidden:
            if f.lower() in txt.lower(): hits.append(f'{p.name}:{f}')
ck('common_side_firewall',len(hits)==0,' | '.join(hits),'0 forbidden input/output references')
rep=pd.DataFrame(checks,columns=['check_id','status','observed','expected','detail']); write_tsv(rep,OUT/'MAR02_postflight_checks.tsv'); bad=rep[rep.status!='PASS']
if len(bad): print(bad.to_string(index=False)); raise SystemExit(f'HOLD: postflight failed {len(bad)} checks')
# Freeze common-side hypotheses.
def sethash(vals): return hash_gene_set(vals)
primary_units=[]
for uid,g in dm[pd.to_numeric(dm.Expanded_flag)==1].groupby('unit_id'):
    genes=sorted(g.HGNC_id.astype(str)); comp=sorted(dm[(dm.unit_id==uid)&(pd.to_numeric(dm.Compact_flag)==1)].HGNC_id.astype(str)); primary_units.append({'unit_id':str(uid),'Expanded_n':len(genes),'Expanded_hash':sethash(genes),'Compact_n':len(comp),'Compact_hash':sethash(comp)})
standalones=sorted(st.pathway_key.astype(str).tolist())
lock={'stage':'MAR02_COMMON_HYPOTHESIS_LOCK','version':VERSION,'freeze_timestamp_utc':utc(),'association_results_read':'NO','gwas_inputs_used':'NO','magma_gene_window_kb':10,'families':{'A_primary_higher_order_Expanded':primary_units,'B_secondary_standalones':standalones,'C_sensitivity_Compact':[{'unit_id':x['unit_id'],'Compact_n':x['Compact_n'],'Compact_hash':x['Compact_hash']} for x in primary_units]},'input_hashes':json.loads((OUT/'MAR02_preanalysis_lock.json').read_text())['input_hashes'],'canonical_output_hashes':{p.name:sha256_file(p) for p in sorted(OUT.glob('MAR02_raw26_*')) if p.is_file()}}
atomic_json(OUT/'MAR02_COMMON_HYPOTHESIS_LOCK.json',lock)
atomic_text(OUT/'MAR02_session_versions.txt','\n'.join([f'timestamp_utc={utc()}',f'python={sys.version.replace(os.linesep," ")}',f'platform={platform.platform()}',f'pandas={pd.__version__}']))
summary_cmp=read_tsv(OUT/'MAR02_raw26_vs_legacy_mdv5_summary.tsv'); allsame=summary_cmp.loc[summary_cmp.component=='ALL_DOWNSTREAM_MEMBERSHIP','exact_same'].iloc[0]
atomic_text(OUT/'MAR02_validation_summary.txt',f'MAR02 RAW26 COMPRESSION VALIDATION SUMMARY\nstatus=PASS\ntier1_n={n}\nrepresentative_n={len(reps)}\nhigher_order_unit_n={len(units)}\nstandalone_or_unstable_n={len(st)}\nall_downstream_membership_exact_same_vs_legacy={bool(int(allsame))}\nassociation_results_read=NO\ngwas_inputs_used=NO\nnext_stage=MAR03')
# Analysis lock + checksums. PASS is last.
files=[p for p in sorted(OUT.iterdir()) if p.is_file() and p.name not in {'MAR02_checksums.sha256','MAR02_PASS.txt'}]
analysis_lock={'stage':'MAR02','version':VERSION,'freeze_timestamp_utc':utc(),'preanalysis_lock_sha256':sha256_file(OUT/'MAR02_preanalysis_lock.json'),'common_hypothesis_lock_sha256':sha256_file(OUT/'MAR02_COMMON_HYPOTHESIS_LOCK.json'),'outputs':{p.name:sha256_file(p) for p in files}}
atomic_json(OUT/'MAR02_analysis_lock.json',analysis_lock)
files=[p for p in sorted(OUT.iterdir()) if p.is_file() and p.name not in {'MAR02_checksums.sha256','MAR02_PASS.txt'}]
atomic_text(OUT/'MAR02_checksums.sha256','\n'.join(f'{sha256_file(p)}  {p.relative_to(MA)}' for p in files))
atomic_text(OUT/'MAR02_PASS.txt',f'status=PASS\ntechnical_status=PASS\nstage=MAR02_RAW26_COMPRESSION\nversion={VERSION}\ntier1_n={n}\nrepresentative_n={len(reps)}\nhigher_order_unit_n={len(units)}\nstandalone_or_unstable_n={len(st)}\nall_downstream_membership_exact_same_vs_legacy={bool(int(allsame))}\nassociation_results_read=NO\ngwas_inputs_used=NO\ncommon_hypothesis_lock_sha256={sha256_file(OUT/"MAR02_COMMON_HYPOTHESIS_LOCK.json")}\nanalysis_lock_sha256={sha256_file(OUT/"MAR02_analysis_lock.json")}')
print(f'MAR02_POSTFLIGHT=PASS\nMAR02_PASS={OUT/"MAR02_PASS.txt"}\nALL_DOWNSTREAM_EXACT_SAME={bool(int(allsame))}')
