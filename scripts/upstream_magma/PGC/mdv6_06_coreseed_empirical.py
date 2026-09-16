#!/usr/bin/env python3
import argparse, gzip, json
from pathlib import Path
import numpy as np, pandas as pd
from mdv6_common import *

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); a=ap.parse_args()
    cfg=load_cfg(a.config); out=out_dir(cfg); root=Path(cfg["project_root"])
    require_file(out/cfg["outputs"]["primary_gsa_pass"])
    g=pd.read_csv(out/cfg["outputs"]["primary_gene_results"],sep="\t")
    if "ZSTAT" not in g.columns: raise RuntimeError("Primary gene results missing ZSTAT")
    g["ZSTAT"]=pd.to_numeric(g.ZSTAT,errors="coerce")
    zsym=dict(zip(g.HGNC_symbol,g.ZSTAT))
    u=pd.read_csv(root_path(cfg,cfg["upstream"]["gene_universe"]),sep="\t")
    core=u.loc[pd.to_numeric(u.CoreSeed_flag,errors="coerce").fillna(0).eq(1),"HGNC_symbol"].tolist()
    obs=np.array([zsym[x] for x in core if x in zsym and np.isfinite(zsym[x])],float); nobs=len(obs)
    if nobs < int(np.ceil(cfg["expected"]["coreseed_n"]*cfg["mapping_qc"]["coreseed_gene_coverage_min"])): raise RuntimeError("Too few mapped CoreSeed genes for empirical test")
    pseudo=pd.read_csv(root_path(cfg,cfg["upstream"]["mdv4_pseudo_sets"]),sep="\t")
    if len(pseudo)!=cfg["primary_inference"]["coreseed_empirical_B"]: raise RuntimeError("MDV4 pseudo set count drift")
    stats=[]
    for r in pseudo.itertuples(index=False):
        genes=list(r)[1:]
        zz=np.array([zsym[x] for x in genes if x in zsym and np.isfinite(zsym[x])],float)
        stats.append((int(r.set_id),len(zz),float(np.mean(zz)) if len(zz) else np.nan,float(np.median(zz)) if len(zz) else np.nan))
    ns=pd.DataFrame(stats,columns=["set_id","mapped_n","mean_Z","median_Z"])
    ns.to_csv(out/cfg["outputs"]["coreseed_null_stats"],sep="\t",index=False,compression="gzip")
    exact=ns[ns.mapped_n==nobs].copy(); minvalid=cfg["primary_inference"]["coreseed_empirical_min_exact_mapped_null_sets"]
    if len(exact)<minvalid: raise RuntimeError(f"Only {len(exact)} pseudo sets have exact mapped_n={nobs}; minimum={minvalid}. Do not change null after viewing results; audit mapping first.")
    rows=[]
    for stat,obsval in [("mean_Z",float(obs.mean())),("median_Z",float(np.median(obs)))]:
        arr=exact[stat].to_numpy(float); exceed=int(np.sum(arr>=obsval)); p=(1+exceed)/(len(arr)+1); percentile=float(np.mean(arr<=obsval))
        rows.append({"statistic":stat,"CoreSeed_mapped_n":nobs,"observed":obsval,"exact_size_null_sets":len(arr),"null_mean":float(np.mean(arr)),"null_sd":float(np.std(arr,ddof=1)),"exceed_n":exceed,"P_emp_greater":p,"percentile":percentile,"null_source":"frozen MDV4 10k structure-matched pseudo-CoreSeeds"})
    res=pd.DataFrame(rows); res.to_csv(out/cfg["outputs"]["coreseed_empirical"],sep="\t",index=False)
    write_text(out/cfg["outputs"]["coreseed_empirical_pass"],f"status=PASS\nstage=MDV6_CORESEED_EMPIRICAL\nCoreSeed_mapped_n={nobs}\nexact_size_null_sets={len(exact)}\nnull_reused_from=MDV4_FROZEN_10000\n")
    print(res.to_string(index=False))
if __name__=="__main__": main()
