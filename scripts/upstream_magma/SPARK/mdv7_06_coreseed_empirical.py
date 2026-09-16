#!/usr/bin/env python3
import argparse
import numpy as np, pandas as pd
from mdv7_common import *

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); a=ap.parse_args(); cfg=load_cfg(a.config); root=Path(cfg["project_root"]); out=out_dir(cfg); require_file(out/cfg["outputs"]["expanded_pass"])
    g=pd.read_csv(out/cfg["outputs"]["primary_gene_results"],sep="\t"); g["ZSTAT"]=pd.to_numeric(g.ZSTAT,errors="coerce"); zsym=dict(zip(g.HGNC_symbol,g.ZSTAT)); u=pd.read_csv(root_path(cfg,cfg["upstream"]["gene_universe"]),sep="\t"); core=u.loc[pd.to_numeric(u.CoreSeed_flag,errors="coerce").fillna(0).eq(1),"HGNC_symbol"].tolist(); obs=np.array([zsym[x] for x in core if x in zsym and np.isfinite(zsym[x])],float); nobs=len(obs)
    if nobs<int(np.ceil(cfg["expected"]["coreseed_n"]*cfg["mapping_qc"]["coreseed_gene_coverage_min"])): raise RuntimeError("Too few mapped CoreSeed genes")
    pseudo=pd.read_csv(root_path(cfg,cfg["upstream"]["mdv4_pseudo_sets"]),sep="\t");
    if len(pseudo)!=cfg["coreseed_inference"]["empirical_B"]: raise RuntimeError("Frozen pseudo-set count drift")
    rows=[]
    for r in pseudo.itertuples(index=False):
        zz=np.array([zsym[x] for x in list(r)[1:] if x in zsym and np.isfinite(zsym[x])],float); rows.append((int(r.set_id),len(zz),float(np.mean(zz)) if len(zz) else np.nan,float(np.median(zz)) if len(zz) else np.nan))
    ns=pd.DataFrame(rows,columns=["set_id","mapped_n","mean_Z","median_Z"]); ns.to_csv(out/cfg["outputs"]["coreseed_null_stats"],sep="\t",index=False,compression="gzip"); exact=ns[ns.mapped_n==nobs];
    if len(exact)<cfg["coreseed_inference"]["empirical_min_exact_mapped_null_sets"]: raise RuntimeError("Insufficient exact-size frozen pseudo-CoreSeeds")
    rr=[]
    for stat,val in [("mean_Z",float(obs.mean())),("median_Z",float(np.median(obs)))]:
        arr=exact[stat].to_numpy(float); exceed=int(np.sum(arr>=val)); rr.append({"statistic":stat,"CoreSeed_mapped_n":nobs,"observed":val,"exact_size_null_sets":len(arr),"null_mean":float(arr.mean()),"null_sd":float(arr.std(ddof=1)),"exceed_n":exceed,"P_emp_greater":(1+exceed)/(len(arr)+1),"percentile":float(np.mean(arr<=val)),"null_source":"frozen MDV4 10k structure-matched pseudo-CoreSeeds"})
    pd.DataFrame(rr).to_csv(out/cfg["outputs"]["coreseed_empirical"],sep="\t",index=False); write_text(out/cfg["outputs"]["coreseed_empirical_pass"],f"status=PASS\nstage=MDV7_CORESEED_EMPIRICAL\nCoreSeed_mapped_n={nobs}\nexact_size_null_sets={len(exact)}\nnull_reused_from=MDV4_FROZEN_10000\nrole=SUPPORTIVE_NOT_PRIMARY\n")
    print("MDV7_CORESEED_EMPIRICAL_PASS")
if __name__=="__main__": main()
