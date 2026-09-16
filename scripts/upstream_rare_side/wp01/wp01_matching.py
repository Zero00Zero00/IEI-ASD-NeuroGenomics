#!/usr/bin/env python3
from pathlib import Path
import sys,json,hashlib,numpy as np,pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from wp01_common import *
from wp01_lineage import verify_stable

FEATURES=["z_log_gene_length","z_GC","z_log_annotation"]

def verify_alignment(out):
    req=[out/"WP01_ALIGNMENT_PASS.txt",out/"WP01_ALIGNMENT_LOCK.json",out/"WP01_model_frame.tsv.gz"]
    if not all(p.exists() for p in req): raise RuntimeError("WP01 alignment checkpoint missing")
    lock=json.loads((out/"WP01_ALIGNMENT_LOCK.json").read_text())
    if lock.get("status")!="PASS_CANDIDATE": raise RuntimeError("Alignment lock not PASS_CANDIDATE")
    if sha256_file(out/"WP01_model_frame.tsv.gz")!=lock["model_frame_sha256"]:
        raise RuntimeError("Canonical model-frame hash changed after alignment")
    return lock

def verify_stage1_artifacts(out):
    obj=json.loads((out/"WP01_stage1_lock.json").read_text())
    for f,h in obj["files"].items():
        p=out/f
        if not p.exists() or sha256_file(p)!=h:
            raise RuntimeError(f"Stage1 frozen artifact changed after STOP-B lock: {f}")
    return obj

def build_pools(df,K):
    anchors=df[df.raw26_flag==1].copy()
    controls=df[df.raw26_flag==0].copy()
    pools=[]
    for _,a in anchors.sort_values("gene").iterrows():
        q=int(a.length_quintile)
        cand=controls[controls.length_quintile==q].copy()
        cand["distance"]=np.sqrt(sum((cand[c]-a[c])**2 for c in FEATURES))
        cand=cand.sort_values(["distance","gene"]).head(K)
        if len(cand)<K:
            remaining=controls[~controls.gene.isin(cand.gene)].copy()
            remaining["qdist"]=(remaining.length_quintile-q).abs()
            remaining["distance"]=np.sqrt(sum((remaining[c]-a[c])**2 for c in FEATURES))
            remaining=remaining.sort_values(["qdist","distance","gene"]).head(K-len(cand))
            cand=pd.concat([cand,remaining],ignore_index=True)
        if len(cand)<K: raise RuntimeError(f"{a.gene}: only {len(cand)} candidates for K={K}")
        for rank,(_,r) in enumerate(cand.iterrows(),1):
            pools.append({"anchor_gene":a.gene,"candidate_gene":r.gene,"rank":rank,
                          "distance":r.distance,"anchor_length_quintile":q,
                          "candidate_length_quintile":int(r.length_quintile)})
    return pd.DataFrame(pools)

def generate_sets(df,pools,n_sets,seed):
    rng=np.random.default_rng(seed)
    anchors=sorted(pools.anchor_gene.unique())
    by={a:pools[pools.anchor_gene==a].sort_values("rank").candidate_gene.to_numpy() for a in anchors}
    records=[]; seen=set(); set_rows=[]; attempts=0
    raw=df[df.raw26_flag==1]
    idx=df.set_index("gene")
    bal=[]
    while len(seen)<n_sets:
        attempts+=1
        if attempts>n_sets*100: raise RuntimeError("Too many attempts generating unique matched sets")
        used=set(); assign={}; ok=True
        for a in rng.permutation(anchors):
            avail=np.array([g for g in by[a] if g not in used],dtype=object)
            if len(avail)==0: ok=False; break
            g=rng.choice(avail); assign[a]=g; used.add(g)
        if not ok: continue
        h=exact_set_hash(used)
        if h in seen: continue
        seen.add(h); rid=len(seen)
        set_rows.append({"replicate_id":rid,"set_sha256":h,"gene_set":";".join(sorted(used))})
        for a in anchors: records.append({"replicate_id":rid,"anchor_gene":a,"control_gene":assign[a]})
        controls=idx.loc[list(used)]
        for col in FEATURES:
            bal.append({"replicate_id":rid,"covariate":col,"SMD":smd_pooled(controls[col],raw[col])})
    return pd.DataFrame(records),pd.DataFrame(set_rows),pd.DataFrame(bal),attempts

def balance_summary(bal,c,K,n,attempts):
    z=bal.assign(absSMD=bal.SMD.abs()).groupby("covariate").absSMD.agg(
        median_abs_SMD="median",p95_abs_SMD=lambda x:x.quantile(.95),max_abs_SMD="max").reset_index()
    z["median_gate"]=float(c["rules"]["matched_balance_median_gate"])
    z["p95_gate"]=float(c["rules"]["matched_balance_p95_gate"])
    z["status"]=np.where((z.median_abs_SMD<=z.median_gate)&(z.p95_abs_SMD<=z.p95_gate),"PASS","FAIL")
    z["attempts"]=attempts; z["accepted_sets"]=n; z["K"]=K
    return z

def main():
    if len(sys.argv)<3: raise SystemExit("usage: wp01_matching.py CH3_ROOT pilot|full")
    root=Path(sys.argv[1]).resolve(); mode=sys.argv[2]
    ma=root/"12_molecular_autism_revision"; out=ma/"01_anchor_ingest/WP01_raw26_primary_v1"
    locksha=sha256_file(out/"WP01_preanalysis_lock.json")
    require_release(out/"WP01_STOPA_RELEASE.txt",locksha,"A")
    verify_stable(ma,out,locksha)
    align=verify_alignment(out)
    c=read_contract(ma); seed=int(c["seeds"]["matched_sets"])
    df=pd.read_csv(out/"WP01_model_frame.tsv.gz",sep="\t",compression="gzip")
    required=["gene","raw26_flag","length_quintile"]+FEATURES
    miss=[x for x in required if x not in df.columns]
    if miss: raise RuntimeError("Canonical matching columns missing: "+",".join(miss))
    for col in FEATURES+["length_quintile","raw26_flag"]:
        df[col]=pd.to_numeric(df[col],errors="raise")
    df["raw26_flag"]=df.raw26_flag.astype(int)
    if len(df)!=19267 or df.gene.nunique()!=19267 or int(df.raw26_flag.sum())!=26:
        raise RuntimeError("Canonical model-frame identity/count failure")

    if mode=="pilot":
        K=int(c["rules"]["matched_primary_K"]); n=int(c["rules"]["matched_pilot_B"])
        pools=build_pools(df,K)
        pools.to_csv(out/"WP01_raw26_candidate_pools_K50.tsv.gz",sep="\t",index=False,compression="gzip")
        ass,sets,bal,attempts=generate_sets(df,pools,n,seed)
        ass.to_csv(out/"WP01_raw26_pilot_assignments_1k.tsv.gz",sep="\t",index=False,compression="gzip")
        sets.to_csv(out/"WP01_raw26_pilot_sets_1k.tsv.gz",sep="\t",index=False,compression="gzip")
        bal.to_csv(out/"WP01_raw26_pilot_balance_1k.tsv.gz",sep="\t",index=False,compression="gzip")
        summary=balance_summary(bal,c,K,n,attempts)
        summary.to_csv(out/"WP01_raw26_pilot_balance_summary.tsv",sep="\t",index=False)
        if len(sets)!=1000 or sets.set_sha256.nunique()!=1000: raise RuntimeError("K50 pilot sets not 1000 unique")
        if (summary.status!="PASS").any(): raise RuntimeError("pilot primary matching balance gate failed")
        print(summary.to_string(index=False))
        return

    if mode!="full": raise RuntimeError("mode must be pilot or full")
    stage=verify_stage1_artifacts(out)
    stage1sha=sha256_file(out/"WP01_stage1_lock.json")
    require_release(out/"WP01_STOPB_RELEASE.txt",locksha,"B","stage1_lock_sha256",stage1sha)

    # Primary K50 / B=10,000. It must reproduce the pilot as an exact prefix.
    K=int(c["rules"]["matched_primary_K"]); n=int(c["rules"]["matched_full_B"])
    pools=pd.read_csv(out/"WP01_raw26_candidate_pools_K50.tsv.gz",sep="\t",compression="gzip")
    ass,sets,bal,attempts=generate_sets(df,pools,n,seed)
    pilot=pd.read_csv(out/"WP01_raw26_pilot_sets_1k.tsv.gz",sep="\t",compression="gzip").sort_values("replicate_id")
    prefix=sets.sort_values("replicate_id").head(1000)
    if not np.array_equal(pilot.set_sha256.to_numpy(),prefix.set_sha256.to_numpy()):
        raise RuntimeError("K50 full run does not reproduce the frozen 1,000-set pilot prefix")
    ass.to_csv(out/"WP01_raw26_full_assignments_10k.tsv.gz",sep="\t",index=False,compression="gzip")
    sets.to_csv(out/"WP01_raw26_full_sets_10k.tsv.gz",sep="\t",index=False,compression="gzip")
    bal.to_csv(out/"WP01_raw26_full_balance_10k.tsv.gz",sep="\t",index=False,compression="gzip")
    s50=balance_summary(bal,c,K,n,attempts)
    s50.to_csv(out/"WP01_raw26_full_balance_summary.tsv",sep="\t",index=False)
    if len(sets)!=10000 or sets.set_sha256.nunique()!=10000 or (s50.status!="PASS").any():
        raise RuntimeError("K50 full matching QC failed")

    # Prespecified K20 sensitivity. Use the SAME frozen matching seed; only K changes.
    K20=int(c["rules"]["matched_sensitivity_K"])
    p20=build_pools(df,K20)
    p20.to_csv(out/"WP01_raw26_candidate_pools_K20.tsv.gz",sep="\t",index=False,compression="gzip")
    a20,s20,b20,att20=generate_sets(df,p20,n,seed)
    a20.to_csv(out/"WP01_raw26_K20_assignments_10k.tsv.gz",sep="\t",index=False,compression="gzip")
    s20.to_csv(out/"WP01_raw26_K20_sets_10k.tsv.gz",sep="\t",index=False,compression="gzip")
    b20.to_csv(out/"WP01_raw26_K20_balance_10k.tsv.gz",sep="\t",index=False,compression="gzip")
    s20sum=balance_summary(b20,c,K20,n,att20)
    s20sum.to_csv(out/"WP01_raw26_K20_balance_summary.tsv",sep="\t",index=False)
    if len(s20)!=10000 or s20.set_sha256.nunique()!=10000 or (s20sum.status!="PASS").any():
        raise RuntimeError("K20 sensitivity matching QC failed")

    print("K50 primary:")
    print(s50.to_string(index=False))
    print("\nK20 sensitivity:")
    print(s20sum.to_string(index=False))

if __name__=="__main__":
    main()
