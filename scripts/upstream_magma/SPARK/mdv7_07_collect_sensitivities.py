#!/usr/bin/env python3
import argparse
import pandas as pd
from mdv7_common import *

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); a=ap.parse_args(); cfg=load_cfg(a.config); out=out_dir(cfg); require_file(out/cfg["outputs"]["compact_pass"])
    scenarios=[("compact","secondary_compact"),("w0","sensitivity_only"),("w50","sensitivity_only"),("mhc","sensitivity_only"),("remove_coreseed","sensitivity_only"),("alltier1","sensitivity_only"),("nested","sensitivity_only")]; rows=[]
    for sc,role in scenarios:
        sd=out/"gsa"/sc; mf=sd/"marginal.tsv"; cf=sd/"conditional.tsv"; require_file(mf); m=pd.read_csv(mf,sep="\t")
        for _,r in m.iterrows(): rows.append({"scenario":sc,"model":"marginal","unit_id":r.unit_id,"NGENES":r.NGENES,"BETA":r.BETA,"SE":r.SE,"P":r.P,"q_BH":r.q_BH,"role":role})
        if cf.exists():
            c=pd.read_csv(cf,sep="\t")
            for _,r in c.iterrows(): rows.append({"scenario":sc,"model":"conditional","unit_id":r.unit_id,"NGENES":r.NGENES,"BETA":r.BETA,"SE":r.SE,"P":r.P,"q_BH":r.q_BH,"role":role})
    z=pd.DataFrame(rows); prim=pd.read_csv(out/cfg["outputs"]["expanded_summary"],sep="\t"); pm=dict(zip(prim.unit_id,prim.beta_marginal)); pc=dict(zip(prim.unit_id,prim.beta_conditional))
    def sdir(r):
        ref=pm.get(r.unit_id) if r.model=="marginal" else pc.get(r.unit_id)
        if ref is None or pd.isna(ref): return pd.NA
        return (r.BETA==0 and ref==0) or (r.BETA*ref>0)
    z["same_direction_as_expanded_primary"]=z.apply(sdir,axis=1); z.to_csv(out/cfg["outputs"]["sensitivity_summary"],sep="\t",index=False); write_text(out/cfg["outputs"]["sensitivity_pass"],"status=PASS\nstage=MDV7_SENSITIVITY\ncompact_role=SECONDARY_7_SET_FAMILY\nother_scenarios=ROBUSTNESS_ONLY_NEVER_RESCUE_EXPANDED_PRIMARY\n")
    print("MDV7_SENSITIVITY_PASS")
if __name__=="__main__": main()
