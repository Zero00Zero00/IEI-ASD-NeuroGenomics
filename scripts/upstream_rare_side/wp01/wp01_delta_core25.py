#!/usr/bin/env python3
from pathlib import Path
import sys,numpy as np,pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from wp01_common import *
from wp01_lineage import verify_stable

def main():
    root=Path(sys.argv[1]).resolve(); ma=root/"12_molecular_autism_revision"
    out=ma/"01_anchor_ingest/WP01_raw26_primary_v1"
    verify_stable(ma,out,sha256_file(out/"WP01_preanalysis_lock.json"))
    raw=pd.read_csv(out/"WP01_raw26_tier.tsv.gz",sep="\t")
    core=pd.read_csv(root/"04_bias_empirical/MDV4_empirical/MDV4_tier_all.tsv.gz",sep="\t",compression="gzip")
    keep=["source","pathway_id","OR_Firth","q_Firth","P_emp","q_emp","firth_pass","emp_pass","Tier1_flag"]
    miss=[x for x in keep if x not in core.columns]
    if miss: raise RuntimeError("Legacy MDV4 binary-gate columns missing: "+",".join(miss))
    c=core[keep].copy().rename(columns={
      "OR_Firth":"core25_OR_Firth","q_Firth":"core25_q_Firth","P_emp":"core25_P_emp",
      "q_emp":"core25_q_emp","firth_pass":"core25_firth_pass",
      "emp_pass":"core25_emp_pass","Tier1_flag":"core25_Tier1_flag"
    })
    m=raw.merge(c,on=["source","pathway_id"],how="left",validate="one_to_one")
    for col in ["core25_firth_pass","core25_emp_pass","core25_Tier1_flag"]:
        m[col]=m[col].map(truthy)
    m["raw26_firth_pass"]=(m.OR_Firth>1)&(m.q_Firth<0.05)
    m["raw26_emp_pass"]=m.q_emp<0.05
    m["raw26_Tier1_flag"]=m.raw26_firth_pass&m.raw26_emp_pass
    m["delta_logOR"]=np.log(pd.to_numeric(m.OR_Firth,errors="coerce"))-np.log(pd.to_numeric(m.core25_OR_Firth,errors="coerce"))
    m["delta_q_Firth"]=pd.to_numeric(m.q_Firth,errors="coerce")-pd.to_numeric(m.core25_q_Firth,errors="coerce")
    m["delta_q_emp"]=pd.to_numeric(m.q_emp,errors="coerce")-pd.to_numeric(m.core25_q_emp,errors="coerce")
    m["Firth_gate_changed"]=m.raw26_firth_pass!=m.core25_firth_pass
    m["Emp_gate_changed"]=m.raw26_emp_pass!=m.core25_emp_pass
    m["Tier1_gate_changed"]=m.raw26_Tier1_flag!=m.core25_Tier1_flag
    m.to_csv(out/"WP01_raw26_vs_core25.tsv.gz",sep="\t",index=False,compression="gzip")

    s=pd.DataFrame([{
      "raw26_Firth_positive_n":int(m.raw26_firth_pass.sum()),
      "core25_Firth_positive_n":int(m.core25_firth_pass.sum()),
      "raw26_emp_positive_n":int(m.raw26_emp_pass.sum()),
      "core25_emp_positive_n":int(m.core25_emp_pass.sum()),
      "raw26_Tier1_n":int(m.raw26_Tier1_flag.sum()),
      "core25_Tier1_n":int(m.core25_Tier1_flag.sum()),
      "Firth_gate_changed_n":int(m.Firth_gate_changed.sum()),
      "Emp_gate_changed_n":int(m.Emp_gate_changed.sum()),
      "Tier1_gate_changed_n":int(m.Tier1_gate_changed.sum()),
      "max_abs_delta_logOR":float(np.nanmax(np.abs(m.delta_logOR))),
      "max_abs_delta_q_Firth":float(np.nanmax(np.abs(m.delta_q_Firth))),
      "max_abs_delta_q_emp":float(np.nanmax(np.abs(m.delta_q_emp)))
    }])
    s.to_csv(out/"WP01_raw26_vs_core25_summary.tsv",sep="\t",index=False)
    print(s.to_string(index=False))

if __name__=="__main__":
    main()
