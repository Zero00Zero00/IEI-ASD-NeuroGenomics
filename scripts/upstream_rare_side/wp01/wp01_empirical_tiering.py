#!/usr/bin/env python3
from pathlib import Path
import sys,math,numpy as np,pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from wp01_common import *
from wp01_lineage import verify_stable

def empirical_for_sets(root,out,sets_file,label):
    mem=load_pathway_membership(root)
    meta=mem[["source","pathway_id","pathway_name"]].drop_duplicates(["source","pathway_id"]).reset_index(drop=True)
    meta["path_index"]=np.arange(len(meta))
    key_to_idx={(r.source,r.pathway_id):int(r.path_index) for r in meta.itertuples()}
    gene_paths={}
    for (g),grp in mem.groupby("gene"):
        gene_paths[g]=np.array([key_to_idx[(s,p)] for s,p in zip(grp.source,grp.pathway_id)],dtype=np.int32)
    path_size=mem.groupby(["source","pathway_id"]).gene.nunique()
    psize=np.array([path_size.loc[(r.source,r.pathway_id)] for r in meta.itertuples()],dtype=int)
    raw=pd.read_csv(out/"WP01_raw26_gene_set.tsv",sep="\t")
    rawgenes=raw.gene.astype(str).tolist(); nset=len(rawgenes)
    obs=np.zeros(len(meta),dtype=np.int16)
    for g in rawgenes:
        if g in gene_paths: obs[gene_paths[g]]+=1
    sets=pd.read_csv(sets_file,sep="\t",compression="gzip")
    exceed=np.zeros(len(meta),dtype=np.int32)
    sumcount=np.zeros(len(meta),dtype=np.int64); maxcount=np.zeros(len(meta),dtype=np.int16)
    for i,r in enumerate(sets.itertuples(index=False),1):
        counts=np.zeros(len(meta),dtype=np.int16)
        for g in str(r.gene_set).split(";"):
            arr=gene_paths.get(g)
            if arr is not None: counts[arr]+=1
        # With fixed universe size, set size and pathway size, the Haldane log-OR is monotone in pathway count.
        exceed += (counts>=obs)
        sumcount += counts
        maxcount=np.maximum(maxcount,counts)
        if i%1000==0: print(f"{label}: {i}/{len(sets)}")
    B=len(sets); pemp=(1+exceed)/(B+1)
    N=19267
    # observed Haldane logOR
    a=obs.astype(float); b=nset-a; c=psize.astype(float)-a; d=N-psize-b
    zero=(a==0)|(b==0)|(c==0)|(d==0)
    aa=a.copy();bb=b.copy();cc=c.copy();dd=d.copy()
    aa[zero]+=0.5;bb[zero]+=0.5;cc[zero]+=0.5;dd[zero]+=0.5
    logOR=np.log((aa*dd)/(bb*cc))
    res=meta[["source","pathway_id","pathway_name"]].copy()
    res["pathway_size"]=psize; res["observed_driver_n"]=obs; res["T_obs_logOR"]=logOR
    res["exceed_n"]=exceed; res["B"]=B; res["P_emp"]=pemp
    res["null_mean_driver_n"]=sumcount/B; res["null_max_driver_n"]=maxcount
    res["q_emp"]=np.nan
    for src,idx in res.groupby("source").groups.items():
        res.loc[idx,"q_emp"]=bh(res.loc[idx,"P_emp"])
    res.to_csv(out/f"WP01_raw26_empirical_{label}.tsv.gz",sep="\t",index=False,compression="gzip")
    return res

def main():
    root=Path(sys.argv[1]).resolve(); ma=root/"12_molecular_autism_revision"; out=ma/"01_anchor_ingest/WP01_raw26_primary_v1"
    verify_stable(ma,out,sha256_file(out/"WP01_preanalysis_lock.json"))
    locksha=sha256_file(out/"WP01_preanalysis_lock.json")
    stage1sha=sha256_file(out/"WP01_stage1_lock.json")
    require_release(out/"WP01_STOPB_RELEASE.txt",locksha,"B","stage1_lock_sha256",stage1sha)
    primary=empirical_for_sets(root,out,out/"WP01_raw26_full_sets_10k.tsv.gz","K50_10k")
    sens1=empirical_for_sets(root,out,out/"WP01_raw26_pilot_sets_1k.tsv.gz","K50_1k")
    sens20=empirical_for_sets(root,out,out/"WP01_raw26_K20_sets_10k.tsv.gz","K20_10k")
    mdv3=pd.read_csv(out/"WP01_raw26_mdv3_all.tsv.gz",sep="\t")
    key=["source","pathway_id"]
    x=mdv3.merge(primary[key+["P_emp","q_emp","observed_driver_n","T_obs_logOR"]],on=key,how="left",validate="one_to_one")
    x=x.merge(sens1[key+["P_emp","q_emp"]].rename(columns={"P_emp":"P_emp_1k","q_emp":"q_emp_1k"}),on=key,how="left")
    x=x.merge(sens20[key+["P_emp","q_emp"]].rename(columns={"P_emp":"P_emp_K20","q_emp":"q_emp_K20"}),on=key,how="left")
    x["tier_class"]=np.where((x.OR_Firth>1)&(x.q_Firth<0.05)&(x.q_emp<0.05),"Tier1",
                      np.where((x.q_Firth<0.05)&(x.q_emp>=0.05),"Firth_only",
                      np.where((x.q_Firth>=0.05)&(x.q_emp<0.05),"Empirical_only","Neither")))
    x.to_csv(out/"WP01_raw26_tier.tsv.gz",sep="\t",index=False,compression="gzip")
    x[x.tier_class=="Tier1"].to_csv(out/"WP01_raw26_tier1.tsv",sep="\t",index=False)
    # MDV2 intersection pathways cross-stage
    mdv2=pd.read_csv(root/"03_intersection/03_intersection_primary_hits.tsv",sep="\t")
    s=choose_col(mdv2.columns,SOURCE_ALIASES); p=choose_col(mdv2.columns,PID_ALIASES)
    if not s or not p: raise RuntimeError("MDV2 primary hits schema unresolved")
    mdv2=mdv2.rename(columns={s:"source",p:"pathway_id"})
    cross=mdv2.merge(x,on=key,how="left")
    cross.to_csv(out/"WP01_raw26_crossstage_MDV2.tsv",sep="\t",index=False)
    print("Raw26 Tier1:",(x.tier_class=="Tier1").sum(),"MDV2 retained Tier1:",(cross.tier_class=="Tier1").sum())
if __name__=="__main__": main()
