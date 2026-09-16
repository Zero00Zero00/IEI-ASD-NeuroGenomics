#!/usr/bin/env python3
import argparse, json
from pathlib import Path
import pandas as pd
from mdv7_common import *

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); ap.add_argument("--scenario",required=True,choices=["w10","w0","w50","mhc"]); a=ap.parse_args(); cfg=load_cfg(a.config); root=Path(cfg["project_root"]); out=out_dir(cfg); res=json.load(open(out/cfg["outputs"]["resolved_inputs"])); require_file(out/cfg["outputs"]["harmonization_pass"])
    magma=res["magma_bin"]; ld=res["ld_prefix"]; gl=out/cfg["outputs"]["gene_locations"]
    if a.scenario=="mhc": window=10; pval=out/cfg["outputs"]["spark_pval_mhc"]; sloc=out/cfg["outputs"]["spark_snploc_mhc"]
    else: window=int(a.scenario[1:]); pval=out/cfg["outputs"]["spark_pval"]; sloc=out/cfg["outputs"]["spark_snploc"]
    mdir=out/"magma"; mdir.mkdir(exist_ok=True); pref=mdir/f"SPARK_EUR_{a.scenario}"
    run([magma,"--annotate",f"window={window},{window}","--snp-loc",sloc,"--gene-loc",gl,"--out",pref],str(pref)+".annotate.wrapper.log")
    annot=str(pref)+".genes.annot"; require_file(annot)
    run([magma,"--bfile",ld,"--pval",pval,"ncol=N","--gene-annot",annot,"--gene-model",cfg["magma"]["primary_gene_model"],"--out",pref],str(pref)+".gene.wrapper.log")
    gout=Path(str(pref)+".genes.out"); graw=Path(str(pref)+".genes.raw"); require_file(gout); require_file(graw)
    if a.scenario=="w10":
        g=read_magma_table(gout); mp=pd.read_csv(out/cfg["outputs"]["gene_id_map"],sep="\t",dtype={"MAGMA_GENE":str})
        if "GENE" not in g.columns: raise RuntimeError("MAGMA genes.out missing GENE")
        g["GENE"]=g.GENE.astype(str)
        if "ZSTAT" not in g.columns and "P" in g.columns:
            from scipy.stats import norm
            pv=pd.to_numeric(g.P,errors="coerce").clip(lower=1e-300,upper=1); g["ZSTAT"]=norm.isf(pv)
        g=g.merge(mp,left_on="GENE",right_on="MAGMA_GENE",how="left"); g.to_csv(out/cfg["outputs"]["primary_gene_results"],sep="\t",index=False)
        mapped=set(g.HGNC_id.dropna());
        def readsets(p):
            d={}
            for line in Path(p).read_text().splitlines():
                z=line.split(); d[z[0]]=set(z[1:])
            return d
        exp=readsets(out/cfg["outputs"]["primary_sets"]); core=readsets(out/cfg["outputs"]["coreseed_set"])["CoreSeed"]; mapped_magma=set(g.GENE.astype(str))
        rows=[{"set":"whole_universe","frozen_n":len(mp),"mapped_n":len(mapped),"coverage":len(mapped)/len(mp)},{"set":"CoreSeed","frozen_n":len(core),"mapped_n":len(core&mapped_magma),"coverage":len(core&mapped_magma)/len(core)}]
        for u in cfg["expected"]["unit_ids"]: rows.append({"set":u,"frozen_n":len(exp[u]),"mapped_n":len(exp[u]&mapped_magma),"coverage":len(exp[u]&mapped_magma)/len(exp[u])})
        q=pd.DataFrame(rows); q.to_csv(out/cfg["outputs"]["gene_mapping_qc"],sep="\t",index=False); fail=[]
        if q.loc[q.set=="whole_universe","coverage"].iloc[0]<cfg["mapping_qc"]["overall_magma_gene_coverage_min"]: fail.append("whole_universe")
        if q.loc[q.set=="CoreSeed","coverage"].iloc[0]<cfg["mapping_qc"]["coreseed_gene_coverage_min"]: fail.append("CoreSeed")
        fail+=q.loc[q.set.isin(cfg["expected"]["unit_ids"])&(q.coverage<cfg["mapping_qc"]["per_unit_gene_coverage_min"]),"set"].tolist()
        if fail and cfg["mapping_qc"]["hard_fail_below_minimum"]: raise RuntimeError("MAGMA mapping gate failed: "+",".join(fail))
        write_text(out/cfg["outputs"]["primary_gene_pass"],f"status=PASS\nstage=MDV7_PRIMARY_MAGMA_GENE\nassociation_rows_read=YES\nwindow_kb=10\ngenes_tested={len(g)}\ncoverage_gate=PASS\nsample_size_mode=ncol=N\n")
    print("MDV7_MAGMA_GENE_DONE",a.scenario)
if __name__=="__main__": main()
