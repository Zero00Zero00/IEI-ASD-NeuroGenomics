#!/usr/bin/env python3
import argparse
import pandas as pd
from scipy.stats import norm
from mdv7_common import *

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); a=ap.parse_args(); cfg=load_cfg(a.config); out=out_dir(cfg); require_file(out/cfg["outputs"]["concordance_pass"])
    pg=pd.read_csv(root_path(cfg,cfg["upstream"]["mdv6_primary_marginal"]),sep="\t")[["unit_id","beta_marginal","SE_marginal","P_marginal","q_marginal"]].rename(columns={"beta_marginal":"PGC_beta_marginal","SE_marginal":"PGC_SE_marginal","P_marginal":"PGC_P_marginal","q_marginal":"PGC_q_marginal"})
    sp=pd.read_csv(out/cfg["outputs"]["expanded_marginal"],sep="\t")[["unit_id","beta_marginal","SE_marginal","P_marginal","q_marginal"]].rename(columns={"beta_marginal":"SPARK_beta_marginal","SE_marginal":"SPARK_SE_marginal","P_marginal":"SPARK_P_marginal","q_marginal":"SPARK_q_marginal"})
    z=pg.merge(sp,on="unit_id")
    alpha=float(cfg["precision"]["conservative_alpha"]); power=float(cfg["precision"]["power"]); mult=norm.ppf(1-alpha)+norm.ppf(power)
    z["SPARK_to_PGC_SE_ratio"]=z.SPARK_SE_marginal/z.PGC_SE_marginal
    z["PGC_MDE_positive_beta_80pct"]=mult*z.PGC_SE_marginal
    z["SPARK_MDE_positive_beta_80pct"]=mult*z.SPARK_SE_marginal
    z["alpha_for_MDE"]=alpha; z["power_for_MDE"]=power; z["MDE_role"]="DESCRIPTIVE_PRECISION_ONLY"
    z.to_csv(out/cfg["outputs"]["precision_audit"],sep="\t",index=False)
    txt=("MDV-7 PRECISION AUDIT\n"+"="*45+"\n"+f"conservative_alpha={alpha}\npower={power}\nMDE_multiplier={mult:.8g}\nmedian_SPARK_to_PGC_SE_ratio={z.SPARK_to_PGC_SE_ratio.median():.8g}\nrole=DESCRIPTIVE_PRECISION_ONLY_NEVER_CHANGES_SIGNIFICANCE\n")
    write_text(out/cfg["outputs"]["precision_summary"],txt); write_text(out/cfg["outputs"]["precision_pass"],"status=PASS\nstage=MDV7_PRECISION_AUDIT\nrole=DESCRIPTIVE_ONLY\n")
    print(txt)
if __name__=="__main__": main()
