#!/usr/bin/env python3
from common import *
import importlib, copy

ensure_dirs(); checks=[]
def ck(cid,ok,obs='',exp='',detail=''):
    checks.append([cid,'PASS' if ok else 'HOLD',str(obs),str(exp),str(detail)])

# Package + upstream authorities.
for name,exp in EXPECTED_WP01_HASHES.items():
    p=WP01/name; ok=p.is_file() and sha256_file(p)==exp; ck('wp01_hash::'+name,ok,sha256_file(p) if p.is_file() else 'MISSING',exp)
for path,exp in EXPECTED_INTERFACE_HASHES.items():
    p=Path(path); ok=p.is_file() and sha256_file(p)==exp; ck('interface_hash::'+p.name,ok,sha256_file(p) if p.is_file() else 'MISSING',exp)
ck('wp01_pass',pass_file(WP01/'WP01_PASS.txt'),WP01/'WP01_PASS.txt','PASS')
ck('legacy_mdv5_pass',pass_file(LEGACY/'MDV5_PASS.txt'),LEGACY/'MDV5_PASS.txt','PASS')

# Re-verify frozen checksum chains rather than trusting manifest filenames alone.
for label,man in [('WP01',WP01/'WP01_checksums.sha256'),('LEGACY_MDV5',LEGACY/'MDV5_checksums.sha256')]:
    nfail,fail=verify_checksum_manifest(man,[CH3,MA,man.parent]); ck('checksum_chain::'+label,len(fail)==0,f'verified={nfail};failures={len(fail)}','0 failures',' | '.join(fail[:5]))
# Legacy analysis lock, when present, must still recognize the same frozen parameter config.
lockp=LEGACY/'MDV5_analysis_lock.json'
if lockp.is_file():
    try:
        lo=json.loads(lockp.read_text()); ch=lo.get('code_hashes',{}); expected_cfg=ch.get('config/mdv5_domain_gate_v1_1.yaml'); ck('legacy_config_lock_hash',expected_cfg is None or expected_cfg==sha256_file(LEGACY_CFG),sha256_file(LEGACY_CFG),expected_cfg or 'not-recorded')
    except Exception as e: ck('legacy_analysis_lock_parse',False,e,'readable JSON')
for p in [UNIVERSE,MEMBERSHIP,LEGACY_CFG]: ck('required::'+p.name,p.is_file() and p.stat().st_size>0,p,'non-empty')
# Dependencies.
vers={}; dep_ok=True; dep_err=''
for name in ['pandas','numpy','scipy','igraph','yaml']:
    try:
        m=importlib.import_module(name); vers[name]=getattr(m,'__version__','unknown')
    except Exception as e: dep_ok=False; dep_err+=f'{name}:{e};'
ck('dependencies',dep_ok,json.dumps(vers,sort_keys=True),'all import',dep_err)

# Validate actual legacy config and inherit scientific parameters rather than guessing them.
try:
    base=yaml.safe_load(LEGACY_CFG.read_text())
    ck('rule::redundancy_0.70',abs(float(base['redundancy']['high_redundancy_jaccard'])-0.70)<1e-12,base['redundancy']['high_redundancy_jaccard'],0.70)
    ck('rule::primary_edge_0.20',abs(float(base['clustering']['primary_edge_threshold'])-0.20)<1e-12,base['clustering']['primary_edge_threshold'],0.20)
    sens=sorted(float(x) for x in base['clustering'].get('sensitivity_edge_thresholds',[])); ck('rule::sensitivity_edges',sens==[0.15,0.25,0.30],sens,[0.15,0.25,0.30])
    ck('rule::leiden_runs_100',int(base['clustering']['n_runs'])==100,base['clustering']['n_runs'],100)
    ck('rule::stable_min3',int(base['clustering']['stable_domain']['min_representative_pathways'])==3,base['clustering']['stable_domain']['min_representative_pathways'],3)
    ck('rule::stable_coclust0.70',abs(float(base['clustering']['stable_domain']['median_within_domain_coclustering_min'])-0.70)<1e-12,base['clustering']['stable_domain']['median_within_domain_coclustering_min'],0.70)
    ck('rule::small_exact2',int(base['clustering']['small_module']['exact_representative_pathways'])==2,base['clustering']['small_module']['exact_representative_pathways'],2)
    ck('rule::small_coclust0.70',abs(float(base['clustering']['small_module']['pair_coclustering_min'])-0.70)<1e-12,base['clustering']['small_module']['pair_coclustering_min'],0.70)
except Exception as e:
    base={}; ck('legacy_config_load',False,e,'readable config')

# Raw26 genes and frozen tier1.
try:
    g=read_tsv(WP01/'WP01_raw26_gene_set.tsv'); cg=pick(g,['gene','HGNC_symbol','approved_symbol']); genes=g[cg].astype(str).str.strip();
    ck('raw26_gene_n',len(g)==26,len(g),26); ck('raw26_gene_unique',genes.nunique()==26,genes.nunique(),26)
    for fl in ['IUIS_flag','SFARI_R0_highconf_flag']:
        if fl in g.columns: ck('raw26_flag::'+fl,g[fl].astype(str).str.lower().isin(['true','1','yes']).all(),g[fl].value_counts().to_dict(),'all TRUE')
except Exception as e: g=None; genes=[]; ck('raw26_load',False,e,'26-gene table')
try:
    t=canon_pathways(read_tsv(WP01/'WP01_raw26_tier1.tsv'))
    for c in ['OR_Firth','q_Firth','q_emp']: t[c]=pd.to_numeric(t[c],errors='coerce')
    rule=(t['OR_Firth']>1)&(t['q_Firth']<0.05)&(t['q_emp']<0.05)
    ck('tier1_nonempty',len(t)>0,len(t),'>0'); ck('tier1_unique',t['pathway_key'].nunique()==len(t),t['pathway_key'].nunique(),len(t)); ck('tier1_rule',bool(rule.all()),int(rule.sum()),len(t))
    if 'tier_class' in t.columns: ck('tier1_label',t['tier_class'].astype(str).eq('Tier1').all(),t['tier_class'].value_counts().to_dict(),'Tier1 only')
except Exception as e: t=None; ck('tier1_load',False,e,'readable tier1')

try:
    u0=read_tsv(UNIVERSE); u=canon_universe(u0); m=canon_membership(read_tsv(MEMBERSHIP));
    # frozen structural checks against legacy config expected values
    ex=base.get('expected',{})
    ck('universe_n',len(u)==int(ex.get('universe_n',len(u))),len(u),ex.get('universe_n',len(u)))
    term=m.groupby(['source','pathway_id','pathway_key']).size().reset_index(name='n'); sc=term.groupby('source')['pathway_key'].nunique().to_dict()
    ck('pathway_total',term['pathway_key'].nunique()==int(ex.get('total_pathways_n',term['pathway_key'].nunique())),term['pathway_key'].nunique(),ex.get('total_pathways_n'))
    ck('pathway_go',int(sc.get('GO_BP',0))==int(ex.get('go_terms_n',sc.get('GO_BP',0))),sc.get('GO_BP',0),ex.get('go_terms_n'))
    ck('pathway_reactome',int(sc.get('Reactome',0))==int(ex.get('reactome_terms_n',sc.get('Reactome',0))),sc.get('Reactome',0),ex.get('reactome_terms_n'))
    if g is not None:
        missing=sorted(set(genes)-set(u['HGNC_symbol'])); ck('raw26_all_in_universe',len(missing)==0,len(missing),0,';'.join(missing))
    if t is not None:
        missingp=sorted(set(t['pathway_key'])-set(m['pathway_key'])); ck('tier1_all_in_membership',len(missingp)==0,len(missingp),0,';'.join(missingp[:20]))
except Exception as e: u0=u=m=None; ck('frozen_structures_load',False,e,'readable universe/membership')

# Fail before generating adapters if any authority/scientific rule failed.
rep=pd.DataFrame(checks,columns=['check_id','status','observed','expected','detail']); write_tsv(rep,OUT/'MAR02_preflight_report.tsv')
bad=rep[rep.status!='PASS']
if len(bad):
    print(bad.to_string(index=False)); raise SystemExit(f'HOLD: MAR02 preflight failed {len(bad)} checks; see {OUT}/MAR02_preflight_report.tsv')

# Build raw26 adapter universe. Internal legacy column name CoreSeed_flag is deliberately reused only inside work/;
# canonical MAR02 outputs rename it to raw26_anchor_flag.
rawset=set(genes)
uadapt=u0.copy(); symcol=pick(uadapt,['HGNC_symbol','approved_symbol','symbol','gene']); corecol=pick(uadapt,['CoreSeed_flag','coreseed_flag','CoreSeed']); uadapt[corecol]=uadapt[symcol].astype(str).str.strip().isin(rawset)
write_tsv(uadapt,WORK/'MAR02_raw26_universe_adapter.tsv')
shutil_tier=WP01/'WP01_raw26_tier1.tsv'; (WORK/'MAR02_raw26_tier1_adapter.tsv').write_bytes(shutil_tier.read_bytes())

cfg=copy.deepcopy(base); cfg['stage']='MAR02_RAW26_COMPRESSION'; cfg['version']=VERSION; cfg['project_root']=str(CH3); cfg['output_dir']=str(WORK)
cfg.setdefault('upstream',{}); cfg['upstream']['universe']=str(WORK/'MAR02_raw26_universe_adapter.tsv'); cfg['upstream']['pathway_membership']=str(MEMBERSHIP); cfg['upstream']['mdv4_tier1']=str(WORK/'MAR02_raw26_tier1_adapter.tsv')
cfg.setdefault('expected',{}); cfg['expected']['coreseed_n']=26; cfg['expected']['tier1_n']=len(t); cnt=t.groupby('source')['pathway_key'].nunique().to_dict(); cfg['expected']['tier1_go_n']=int(cnt.get('GO_BP',0)); cfg['expected']['tier1_reactome_n']=int(cnt.get('Reactome',0))
# internal output contract for vendor stages 02-05
cfg['outputs'].update({
'preflight_pass':'MAR02_PREFLIGHT_PASS.txt','tier1_annotated':'internal_tier1_annotated.tsv','jaccard_matrix':'internal_jaccard_matrix.tsv.gz','jaccard_pairs':'internal_jaccard_pairs.tsv.gz','jaccard_pass':'MAR02_JACCARD_PASS.txt',
'redundancy_map':'internal_redundancy_map.tsv','representatives':'internal_representatives.tsv','redundancy_summary':'internal_redundancy_summary.txt','redundancy_pass':'MAR02_REDUNDANCY_PASS.txt',
'network_edges':'internal_network_edges.tsv','network_edges_sensitivity':'internal_network_edges_sensitivity.tsv.gz','leiden_runs':'internal_leiden_runs.tsv.gz','coclustering_matrix':'internal_coclustering_matrix.tsv.gz','coclustering_pairs':'internal_coclustering_pairs.tsv.gz','medoid_runs':'internal_medoid_runs.tsv','primary_partition':'internal_primary_partition.tsv','threshold_sensitivity':'internal_threshold_sensitivity.tsv','leiden_summary':'internal_leiden_summary.txt','leiden_pass':'MAR02_LEIDEN_PASS.txt',
'pathway_to_unit':'internal_pathway_to_unit.tsv','domain_membership':'internal_domain_membership.tsv','domain_summary':'internal_domain_summary.tsv','locus_dominance_audit':'internal_locus_dominance.tsv','standalone_pathways':'internal_standalones.tsv','domain_freeze_summary':'internal_domain_freeze_summary.txt','domain_freeze_complete':'MAR02_DOMAIN_FREEZE_PASS.txt'})
(WORK/'MAR02_adapter_config.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False),encoding='utf-8')
atomic_text(WORK/'MAR02_PREFLIGHT_PASS.txt',f'status=PASS\nstage=MAR02_PREFLIGHT\ntier1_n={len(t)}\nraw26_n=26\ngwas_inputs_used=NO')
resolved={'WP01':str(WP01),'raw26_gene_set':str(WP01/'WP01_raw26_gene_set.tsv'),'raw26_tier1':str(WP01/'WP01_raw26_tier1.tsv'),'universe':str(UNIVERSE),'pathway_membership':str(MEMBERSHIP),'legacy_config':str(LEGACY_CFG),'legacy_mdv5':str(LEGACY),'output':str(OUT),'generated_config':str(WORK/'MAR02_adapter_config.yaml')}
atomic_json(OUT/'MAR02_resolved_inputs.json',resolved)
atomic_text(OUT/'MAR02_preflight_summary.txt',f'MAR02 preflight PASS\ntier1_n={len(t)}\nraw26_n=26\nGO_BP_tier1_n={cnt.get("GO_BP",0)}\nReactome_tier1_n={cnt.get("Reactome",0)}\nlegacy_parameter_source={LEGACY_CFG}\nassociation_results_read=NO')
print(f'MAR02_PREFLIGHT=PASS\nTIER1_N={len(t)}\nRAW26_N=26\nCONFIG={WORK/"MAR02_adapter_config.yaml"}')
