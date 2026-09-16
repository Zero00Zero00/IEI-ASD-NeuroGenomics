#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
from mdv6_common import *

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); a=ap.parse_args(); cfg=load_cfg(a.config); out=out_dir(cfg)
    scenarios=["w0","w50","mhc","compact","remove_coreseed","alltier1","nested"]
    rows=[]
    for sc in scenarios:
        sd=out/"gsa"/sc
        mf=sd/"marginal.tsv"; cf=sd/"conditional.tsv"
        require_file(mf)
        m=pd.read_csv(mf,sep="\t")
        for _,r in m.iterrows():
            rows.append({"scenario":sc,"model":"marginal","unit_id":r.unit_id,"NGENES":r.NGENES,"BETA":r.BETA,"SE":r.SE,"P":r.P,"q_BH":r.q_BH,"role":"sensitivity_only"})
        if cf.exists():
            c=pd.read_csv(cf,sep="\t")
            for _,r in c.iterrows(): rows.append({"scenario":sc,"model":"conditional","unit_id":r.unit_id,"NGENES":r.NGENES,"BETA":r.BETA,"SE":r.SE,"P":r.P,"q_BH":r.q_BH,"role":"sensitivity_only"})
    res=pd.DataFrame(rows)
    # compare same units to primary; nested M04_EXCLUSIVE kept separate.
    prim=pd.read_csv(out/cfg["outputs"]["primary_gsa_summary"],sep="\t")
    pm=dict(zip(prim.unit_id,prim.beta_marginal)); pc=dict(zip(prim.unit_id,prim.beta_conditional))
    def same_dir(r):
        ref=pm.get(r.unit_id) if r.model=="marginal" else pc.get(r.unit_id)
        if ref is None or pd.isna(ref): return pd.NA
        return (r.BETA==0 and ref==0) or (r.BETA*ref>0)
    res["same_direction_as_primary"]=res.apply(same_dir,axis=1)
    res.to_csv(out/cfg["outputs"]["sensitivity_summary"],sep="\t",index=False)
    write_text(out/cfg["outputs"]["sensitivity_pass"],f"status=PASS\nstage=MDV6_SENSITIVITY\nscenarios={';'.join(scenarios)}\nrole=ROBUSTNESS_ONLY_NEVER_RESCUE_PRIMARY\n")
    print("MDV6_SENSITIVITY_PASS")
if __name__=="__main__": main()
