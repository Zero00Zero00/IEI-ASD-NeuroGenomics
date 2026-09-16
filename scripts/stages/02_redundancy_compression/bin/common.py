#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import os, re, csv, json, hashlib, gzip, datetime as dt, subprocess, sys
import pandas as pd
import yaml

CH3=Path(os.environ.get('CH3_ROOT','/home/h3021/chapter3')).resolve()
MA=Path(os.environ.get('MA_ROOT','/home/h3021/chapter3/12_molecular_autism_revision')).resolve()
WP01=MA/'01_anchor_ingest/WP01_raw26_primary_v1'
OUT=MA/'02_raw26_domains/MAR02_raw26_compression_v1'
WORK=OUT/'work'
LOGS=OUT/'logs'
STATE=OUT/'state'
LEGACY=CH3/'05_domains/MDV5_v1_1'
LEGACY_CFG=CH3/'config/mdv5_domain_gate_v1_1.yaml'
UNIVERSE=CH3/'01_gene_universe/01_gene_universe.tsv'
MEMBERSHIP=CH3/'02_pathways/03_pathway_membership.tsv.gz'
PKG=Path(__file__).resolve().parent.parent
VENDOR=PKG/'workflow/legacy_mdv5_snapshot'
VERSION='1.0R1'

EXPECTED_WP01_HASHES={
'WP01_raw26_tier1.tsv':'d0949af3a4af4b7b17810291f16ec83e4d0219259d6258120c3e1460bde70c1d',
'WP01_raw26_tier.tsv.gz':'6c2b6551d434e719bd49f1d8833b3a4807abf77e53917ce27ebddfdbd8512a91',
'WP01_raw26_mdv3_all.tsv.gz':'abc4e7a89864b1fcaeed5831a580837c303a017179c71efe67284a2642911893',
'WP01_raw26_crossstage_MDV2.tsv':'a8c2a27a77b23c5e640264c3dc8f1a684527ee531249721c07b6a0b8a0e8b35a',
'WP01_raw26_gene_set.tsv':'20655a21e2442dd77215cc5010aef56d2e0e293601dfbba975642ead56011198',
'WP01_raw26_vs_core25.tsv.gz':'ad90aa10e889ed3b1e68c10d5145030ce787dc29d2ec6f6657e83965c97330b4',
'WP01_PASS.txt':'c26de83abdc9496cd1301470918ec39ecca1cf389ef0ae1c6504da1800d1c431',
'WP01_checksums.sha256':'8f659357608f9f37c5967dc5925ec5ba113424a41dee8f6e6cef63c8de19cee9',
}
EXPECTED_INTERFACE_HASHES={
'/home/h3021/chapter3/run_mdv5_staged.sh':'98319d9c06189efa29b9382c7d1817d375db56a0ef321e5281115fb3abeacfa8',
'/home/h3021/chapter3/Snakefile.mdv5':'1eb2fe587258d4fce74f7aa51f6c3a253598cc6f0219fa8bcb4e8619c610a6d9',
'/home/h3021/chapter3/05_domains/MDV5_v1_1/MDV5_checksums.sha256':'99526701ba98bbb8aa85f6105560a720fc55a81a4533f9490fdfdecc427e6e64',
}

def ensure_dirs():
    for p in [OUT,WORK,LOGS,STATE]: p.mkdir(parents=True,exist_ok=True)

def sha256_file(p:Path)->str:
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def read_tsv(p:Path, **kwargs):
    return pd.read_csv(p,sep='\t',compression='gzip' if p.suffix=='.gz' else 'infer',low_memory=False,**kwargs)

def write_tsv(df,p:Path,gzip_out=None):
    p.parent.mkdir(parents=True,exist_ok=True)
    if gzip_out is None: gzip_out=p.suffix=='.gz'
    df.to_csv(p,sep='\t',index=False,compression='gzip' if gzip_out else None)

def atomic_text(p:Path,text:str):
    p.parent.mkdir(parents=True,exist_ok=True); tmp=p.with_suffix(p.suffix+'.tmp')
    tmp.write_text(text if text.endswith('\n') else text+'\n',encoding='utf-8'); os.replace(tmp,p)

def atomic_json(p:Path,obj):
    atomic_text(p,json.dumps(obj,indent=2,sort_keys=True,ensure_ascii=False))

def pass_file(p:Path):
    if not p.is_file() or p.stat().st_size==0:return False
    return 'PASS' in p.read_text(encoding='utf-8',errors='replace').upper()

def utc(): return dt.datetime.now(dt.timezone.utc).isoformat()

def norm_source(v):
    s=str(v).strip(); k=re.sub(r'[^a-z0-9]','',s.lower())
    if k in {'go','gobp','geneontologybiologicalprocess','biologicalprocess'}: return 'GO_BP'
    if k in {'reactome','react'}: return 'Reactome'
    return s

def pick(df, aliases, required=True):
    canon={re.sub(r'[^a-z0-9]','',str(c).lower()):c for c in df.columns}
    for a in aliases:
        k=re.sub(r'[^a-z0-9]','',a.lower())
        if k in canon:return canon[k]
    if required: raise RuntimeError(f'Missing column aliases={aliases}; columns={list(df.columns)}')
    return None

def canon_pathways(df):
    df=df.copy(); cs=pick(df,['source','pathway_source','database']); ci=pick(df,['pathway_id','term_id','ID','id']); cn=pick(df,['pathway_name','term_name','name','Pathway'])
    df=df.rename(columns={cs:'source',ci:'pathway_id',cn:'pathway_name'}); df['source']=df['source'].map(norm_source); df['pathway_id']=df['pathway_id'].astype(str).str.strip(); df['pathway_name']=df['pathway_name'].astype(str).str.strip(); df['pathway_key']=df['source']+'::'+df['pathway_id']; return df

def canon_membership(df):
    df=canon_pathways(df); hid=pick(df,['HGNC_id','hgnc_id','HGNC ID']); sym=pick(df,['approved_symbol','HGNC_symbol','symbol','gene']); df=df.rename(columns={hid:'HGNC_id',sym:'approved_symbol'}); df['HGNC_id']=df['HGNC_id'].astype(str).str.strip(); df['approved_symbol']=df['approved_symbol'].astype(str).str.strip(); return df

def canon_universe(df):
    df=df.copy(); hid=pick(df,['HGNC_id','hgnc_id','HGNC ID']); sym=pick(df,['HGNC_symbol','approved_symbol','symbol','gene']); core=pick(df,['CoreSeed_flag','coreseed_flag','CoreSeed']); df=df.rename(columns={hid:'HGNC_id',sym:'HGNC_symbol',core:'CoreSeed_flag'}); df['HGNC_id']=df['HGNC_id'].astype(str).str.strip(); df['HGNC_symbol']=df['HGNC_symbol'].astype(str).str.strip(); return df

def load_generated_config():
    p=WORK/'MAR02_adapter_config.yaml'
    if not p.is_file(): raise SystemExit('HOLD: run preflight first; generated config missing')
    return yaml.safe_load(p.read_text())

def require_release():
    rel=OUT/'MAR02_STOPA_RELEASE.txt'; lock=OUT/'MAR02_preanalysis_lock.json'
    if not pass_file(rel): raise SystemExit('HOLD: manual Gate A not released. Run: bash run_mar02.sh release RELEASE')
    if not lock.is_file(): raise SystemExit('HOLD: preanalysis lock missing')
    txt=rel.read_text(errors='replace'); h=sha256_file(lock)
    if f'preanalysis_lock_sha256={h}' not in txt: raise SystemExit('HOLD: release token does not match current preanalysis lock')

def run_vendor(script:str, config:Path, log:Path):
    env=os.environ.copy(); env['ROOT']=str(CH3); env['PYTHONPATH']=str(VENDOR)+os.pathsep+env.get('PYTHONPATH','')
    cmd=[sys.executable,str(VENDOR/script),'--config',str(config)]
    with open(log,'w',encoding='utf-8') as f:
        cp=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env)
    if cp.returncode!=0:
        tail='\n'.join(log.read_text(encoding='utf-8',errors='replace').splitlines()[-40:])
        raise SystemExit(f'HOLD: {script} failed rc={cp.returncode}. Tail:\n{tail}')

def verify_checksum_manifest(manifest:Path, roots):
    failures=[]; n=0
    if not manifest.is_file(): return 0,[f'MISSING_MANIFEST:{manifest}']
    for line in manifest.read_text(encoding='utf-8',errors='replace').splitlines():
        line=line.strip()
        if not line or line.startswith('#'): continue
        m=re.match(r'^([0-9a-fA-F]{64})\s+[* ]?(.+?)\s*$',line)
        if not m:
            failures.append(f'UNPARSEABLE:{line}'); continue
        n+=1; exp,rel=m.group(1).lower(),m.group(2); rp=Path(rel); cands=[rp] if rp.is_absolute() else [Path(r)/rp for r in roots]+[manifest.parent/rp]
        target=next((x for x in cands if x.is_file()),None)
        if target is None: failures.append(f'MISSING:{rel}'); continue
        act=sha256_file(target)
        if act.lower()!=exp: failures.append(f'MISMATCH:{rel}:{exp}:{act}')
    return n,failures

def hash_gene_set(vals):
    s='\n'.join(sorted(map(str,vals)))+'\n'; return hashlib.sha256(s.encode()).hexdigest()
