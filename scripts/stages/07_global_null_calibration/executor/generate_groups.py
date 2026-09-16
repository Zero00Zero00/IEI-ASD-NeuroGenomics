
from pathlib import Path
import csv, gzip, json, hashlib, os, math, statistics, subprocess, tarfile, shutil, time
from collections import defaultdict, Counter
from datetime import datetime, timezone

CH3 = Path(os.environ.get("CH3_ROOT","/home/h3021/chapter3"))
MA  = Path(os.environ.get("MA_ROOT","/home/h3021/chapter3/12_molecular_autism_revision"))
PKG = Path(__file__).resolve().parent.parent
OUT = MA / "07_chain_calibration" / "MAR07_chain_calibration_v1.0R4p1"
WP01 = MA / "01_anchor_ingest" / "WP01_raw26_primary_v1"
MAR03 = MA / "03_common_retest" / "MAR03_common_variant_retest_v1.0R2"
MAR04 = MA / "04_effect_boundary" / "MAR04_effect_boundary_v1"
MAR05 = MA / "05_magma_collinearity" / "MAR05_collinearity_v1"
MAR06 = MA / "06_matched_null_audit" / "MAR06_matched_null_audit_v1.0R2"
R3 = MA / "07_chain_calibration" / "MAR07_chain_calibration_v1.0R3"

def utc():
    return datetime.now(timezone.utc).isoformat()

def sha256_file(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()

def open_text(p):
    if str(p).endswith(".gz"):
        return gzip.open(p,"rt",encoding="utf-8-sig",errors="replace")
    return open(p,"rt",encoding="utf-8-sig",errors="replace")

def read_tsv(p):
    with open_text(p) as f:
        r=csv.DictReader(f,delimiter="\t")
        return r.fieldnames or [], list(r)

def write_tsv(p,fields,rows):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True)
    with open(p,"w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

def parse_pass(p):
    d={}
    for line in Path(p).read_text(encoding="utf-8",errors="replace").splitlines():
        if "=" in line:
            k,v=line.split("=",1); d[k]=v
    return d

def wilson_ci(k,n,z=1.959963984540054):
    if n<=0: return (float("nan"),float("nan"))
    p=k/n
    den=1+z*z/n
    center=(p+z*z/(2*n))/den
    half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return max(0.0,center-half),min(1.0,center+half)

def ensure_dirs():
    for d in ["state","logs","observed","certification","calibration","source_data","executor_work"]:
        (OUT/d).mkdir(parents=True,exist_ok=True)

import sys, importlib.util, numpy as np, pandas as pd
ensure_dirs()
rep=int(sys.argv[1]); seed=int(sys.argv[2]); outdir=Path(sys.argv[3]); outdir.mkdir(parents=True,exist_ok=True)
resolved=json.loads((OUT/"MAR07_resolved_inputs.json").read_text())
script=Path(resolved["WP01_GLOBAL"])
sys.path.insert(0,str(script.parent))
spec=importlib.util.spec_from_file_location("wp01_global_real",script)
gm=importlib.util.module_from_spec(spec); spec.loader.exec_module(gm)
frame,meta=gm.load_structural_frame(CH3)
q=5
for col in ["log_gene_length","GC","log_annotation"]:
    frame[col+"_bin"]=gm.qbin(frame[col],q)
frame["joint_stratum"]=frame["log_gene_length_bin"].astype(str)+"|"+frame["GC_bin"].astype(str)+"|"+frame["log_annotation_bin"].astype(str)
feat=["z_log_gene_length","z_GC","z_log_annotation"]; raw=["log_gene_length","GC","log_annotation"]
group_sizes={"Overlap_raw":26,"IEI_only":474,"ASD_only":909}
fitted={}
for g in group_sizes:
    mask=(frame.R0_group.astype(str)==g).to_numpy()
    obs,obj,lam,diag,fit=gm.fit_stratified_tilt(frame,mask,"joint_stratum",feat,g)
    fitted[g]=(obs,gm.make_probabilities(lam,obj))

def draw_group(rng,obs,probs,used):
    for tries in range(1,5001):
        picked=[]
        for st,(n,idx,p) in probs.items():
            avail_mask=np.array([i not in used for i in idx],dtype=bool)
            ai=idx[avail_mask]; ap=p[avail_mask]
            if len(ai)<n: break
            ap=ap/ap.sum()
            picked.extend(rng.choice(ai,size=n,replace=False,p=ap).tolist())
        if len(picked)!=len(obs): continue
        P=frame.loc[picked]
        smd={c:gm.smd_pooled(P[c],obs[c]) for c in raw}
        if max(abs(float(x)) for x in smd.values())<=0.20:
            return picked,smd,tries
    raise RuntimeError("Unable to generate balanced mutually exclusive pseudo-group")

ss=np.random.SeedSequence(seed); rngs=[np.random.default_rng(x) for x in ss.spawn(3)]
used=set(); rec=[]; bal=[]
for g,rng in zip(["Overlap_raw","IEI_only","ASD_only"],rngs):
    obs,probs=fitted[g]
    picked,smd,tries=draw_group(rng,obs,probs,used)
    used.update(picked)
    for i in picked: rec.append({"gene":frame.loc[i,"gene"],"R0_group":g})
    for c,v in smd.items(): bal.append({"group":g,"covariate":c,"SMD":v,"draw_attempts":tries})
for i in frame.index:
    if i not in used: rec.append({"gene":frame.loc[i,"gene"],"R0_group":"Background"})
if Counter(z["R0_group"] for z in rec)!=Counter({"Overlap_raw":26,"IEI_only":474,"ASD_only":909,"Background":17858}):
    raise RuntimeError("Pseudo-group counts mismatch")
write_tsv(outdir/"pseudo_groups.tsv",["gene","R0_group"],rec)
write_tsv(outdir/"pseudo_group_balance.tsv",["group","covariate","SMD","draw_attempts"],bal)
print("GROUP_GENERATOR=PASS")
