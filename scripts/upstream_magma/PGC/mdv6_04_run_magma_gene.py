#!/usr/bin/env python3
import argparse, json, os
from pathlib import Path
import pandas as pd
from mdv6_common import *

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); ap.add_argument("--scenario",required=True,choices=["w10","w0","w50","mhc"]); a=ap.parse_args()
    cfg=load_cfg(a.config); out=out_dir(cfg); resolved=json.load(open(out/cfg["outputs"]["resolved_inputs"]))
    require_file(out/cfg["outputs"]["pgc_qc_pass"])
    magma=resolved["magma_bin"]; ld=resolved["ld_prefix"]
    gene_loc=out/cfg["outputs"]["gene_locations"]
    if a.scenario=="mhc": window=10; pval=out/cfg["outputs"]["pgc_pval_mhc"]; sloc=out/cfg["outputs"]["pgc_snploc_mhc"]
    else:
        window=int(a.scenario[1:]); pval=out/cfg["outputs"]["pgc_pval"]; sloc=out/cfg["outputs"]["pgc_snploc"]
    mdir=out/"magma"; mdir.mkdir(exist_ok=True)
    prefix=mdir/f"PGC2019_{a.scenario}"
    # Annotation uses symmetric window in kb.
    run([magma,"--annotate",f"window={window},{window}","--snp-loc",sloc,"--gene-loc",gene_loc,"--out",prefix],str(prefix)+".annotate.wrapper.log")
    annot=str(prefix)+".genes.annot"; require_file(annot)
    run([magma,"--bfile",ld,"--pval",pval,f"N={cfg['pgc2019']['sample_size_primary']}","--gene-annot",annot,"--gene-model",cfg["magma"]["primary_gene_model"],"--out",prefix],str(prefix)+".gene.wrapper.log")
    gout=Path(str(prefix)+".genes.out"); graw=Path(str(prefix)+".genes.raw"); require_file(gout); require_file(graw)
    if a.scenario=="w10":
        g=read_magma_table(gout)
        # enrich with HGNC map
        mp=pd.read_csv(out/cfg["outputs"]["gene_id_map"],sep="\t",dtype={"MAGMA_GENE":str})
        if "GENE" not in g.columns: raise RuntimeError("MAGMA genes.out missing GENE")
        g["GENE"]=g["GENE"].astype(str); z=None
        if "ZSTAT" not in g.columns and "P" in g.columns:
            from scipy.stats import norm
            pv=pd.to_numeric(g["P"],errors="coerce").clip(lower=1e-300,upper=1)
            g["ZSTAT"]=norm.isf(pv)
        g=g.merge(mp,left_on="GENE",right_on="MAGMA_GENE",how="left")
        g.to_csv(out/cfg["outputs"]["primary_gene_results"],sep="\t",index=False)
        universe_n=len(mp); mapped_ids=set(g["HGNC_id"].dropna()); core=set(mp.loc[pd.to_numeric(mp.CoreSeed_flag,errors="coerce").fillna(0).eq(1),"HGNC_id"])
        dm=pd.read_csv(root_path(cfg,cfg["upstream"]["mdv5_domain_membership"]),sep="\t")
        rows=[]
        rows.append({"set":"whole_universe","frozen_n":universe_n,"mapped_n":len(mapped_ids),"coverage":len(mapped_ids)/universe_n})
        rows.append({"set":"CoreSeed","frozen_n":len(core),"mapped_n":len(core & mapped_ids),"coverage":len(core & mapped_ids)/len(core)})
        for uid in cfg["expected"]["unit_ids"]:
            S=set(dm.loc[(dm.domain_id==uid)&(dm.Expanded_flag==1),"HGNC_id"])
            rows.append({"set":uid,"frozen_n":len(S),"mapped_n":len(S & mapped_ids),"coverage":len(S & mapped_ids)/len(S)})
        q=pd.DataFrame(rows); q.to_csv(out/cfg["outputs"]["gene_mapping_qc"],sep="\t",index=False)
        fail=[]
        if q.loc[q.set=="whole_universe","coverage"].iloc[0] < cfg["mapping_qc"]["overall_magma_gene_coverage_min"]: fail.append("whole_universe")
        if q.loc[q.set=="CoreSeed","coverage"].iloc[0] < cfg["mapping_qc"]["coreseed_gene_coverage_min"]: fail.append("CoreSeed")
        fail += q.loc[q.set.isin(cfg["expected"]["unit_ids"]) & (q.coverage<cfg["mapping_qc"]["per_unit_gene_coverage_min"]),"set"].tolist()
        if fail and cfg["mapping_qc"]["hard_fail_below_minimum"]: raise RuntimeError("MAGMA gene mapping coverage gate failed: "+",".join(fail))
        write_text(out/cfg["outputs"]["primary_gene_pass"],f"status=PASS\nstage=MDV6_PRIMARY_MAGMA_GENE\nassociation_results_read=YES\nwindow_kb=10\ngenes_tested={len(g)}\ncoverage_gate=PASS\n")
    print("MAGMA_GENE_DONE",a.scenario)
if __name__=="__main__": main()
