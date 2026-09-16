#!/usr/bin/env python3
import argparse
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from mdv7_common import *

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); a=ap.parse_args(); cfg=load_cfg(a.config); out=out_dir(cfg); require_file(out/cfg["outputs"]["sensitivity_pass"])
    pgm=pd.read_csv(root_path(cfg,cfg["upstream"]["mdv6_primary_marginal"]),sep="\t")[["unit_id","beta_marginal","SE_marginal","P_marginal","q_marginal"]].rename(columns={"beta_marginal":"PGC_beta_marginal","SE_marginal":"PGC_SE_marginal","P_marginal":"PGC_P_marginal","q_marginal":"PGC_q_marginal"})
    pgc=pd.read_csv(root_path(cfg,cfg["upstream"]["mdv6_primary_conditional"]),sep="\t")[["unit_id","beta_conditional","SE_conditional","P_conditional","q_conditional"]].rename(columns={"beta_conditional":"PGC_beta_conditional","SE_conditional":"PGC_SE_conditional","P_conditional":"PGC_P_conditional","q_conditional":"PGC_q_conditional"})
    spm=pd.read_csv(out/cfg["outputs"]["expanded_marginal"],sep="\t")[["unit_id","beta_marginal","SE_marginal","P_marginal","q_marginal"]].rename(columns={"beta_marginal":"SPARK_beta_marginal","SE_marginal":"SPARK_SE_marginal","P_marginal":"SPARK_P_marginal","q_marginal":"SPARK_q_marginal"})
    spc=pd.read_csv(out/cfg["outputs"]["expanded_conditional"],sep="\t")[["unit_id","beta_conditional","SE_conditional","P_conditional","q_conditional"]].rename(columns={"beta_conditional":"SPARK_beta_conditional","SE_conditional":"SPARK_SE_conditional","P_conditional":"SPARK_P_conditional","q_conditional":"SPARK_q_conditional"})
    z=pgm.merge(pgc,on="unit_id").merge(spm,on="unit_id").merge(spc,on="unit_id")
    z["PGC_primary_positive"]=(z.PGC_beta_marginal>0)&(z.PGC_q_marginal<0.05)
    z["SPARK_primary_positive"]=(z.SPARK_beta_marginal>0)&(z.SPARK_q_marginal<0.05)
    z["marginal_sign_concordant"]=(np.sign(z.PGC_beta_marginal)==np.sign(z.SPARK_beta_marginal))
    z["conditional_sign_concordant"]=(np.sign(z.PGC_beta_conditional)==np.sign(z.SPARK_beta_conditional))
    z["cross_cohort_primary_class"]="concordant_primary_null"
    z.loc[(~z.PGC_primary_positive)&z.SPARK_primary_positive,"cross_cohort_primary_class"]="SPARK_only_independent_support_with_heterogeneity"
    z.loc[z.PGC_primary_positive&z.SPARK_primary_positive,"cross_cohort_primary_class"]="both_positive"
    z.loc[z.PGC_primary_positive&(~z.SPARK_primary_positive),"cross_cohort_primary_class"]="PGC_only"
    # Compact family: PGC values came from MDV6 prespecified compact sensitivity; SPARK is MDV7 secondary family.
    pgs=pd.read_csv(root_path(cfg,cfg["upstream"]["mdv6_sensitivity_summary"]),sep="\t")
    pgcomp=pgs[(pgs.scenario=="compact")&(pgs.model=="marginal")][["unit_id","BETA","SE","P","q_BH"]].rename(columns={"BETA":"PGC_compact_beta","SE":"PGC_compact_SE","P":"PGC_compact_P","q_BH":"PGC_compact_q"})
    spcomp=pd.read_csv(out/cfg["outputs"]["compact_marginal"],sep="\t")[["unit_id","beta_marginal","SE_marginal","P_marginal","q_marginal"]].rename(columns={"beta_marginal":"SPARK_compact_beta","SE_marginal":"SPARK_compact_SE","P_marginal":"SPARK_compact_P","q_marginal":"SPARK_compact_q"})
    z=z.merge(pgcomp,on="unit_id",how="left").merge(spcomp,on="unit_id",how="left")
    z["compact_FDR_concordant"]=(z.PGC_compact_beta>0)&(z.PGC_compact_q<0.05)&(z.SPARK_compact_beta>0)&(z.SPARK_compact_q<0.05)
    z.to_csv(out/cfg["outputs"]["cross_cohort"],sep="\t",index=False)
    rho,pv=spearmanr(z.PGC_beta_marginal,z.SPARK_beta_marginal)
    sign_n=int(z.marginal_sign_concordant.sum())
    m02=z[z.unit_id=="M02"].iloc[0]
    if bool(m02.compact_FDR_concordant): m02class="FDR_CONCORDANT_SECONDARY_COMPACT"
    elif float(m02.SPARK_compact_beta)>0 and float(m02.SPARK_compact_P)<0.05: m02class="NOMINAL_DIRECTIONAL_SECONDARY_COMPACT"
    else: m02class="NO_INDEPENDENT_SPARK_SUPPORT_FOR_PGC_COMPACT_SIGNAL"
    txt=("MDV-7 PGC–SPARK CROSS-COHORT CONCORDANCE\n"+"="*62+"\n"+
         f"n_units=7\nmarginal_sign_concordant_n={sign_n}\nSpearman_rho_marginal_beta={rho:.8g}\nSpearman_P_descriptive={pv:.8g}\nSPARK_primary_positive_n={int(z.SPARK_primary_positive.sum())}\nPGC_primary_positive_n={int(z.PGC_primary_positive.sum())}\nM02_compact_class={m02class}\nreplication_language_allowed=NO\nrole=DESCRIPTIVE_CONCORDANCE_AUDIT\n")
    write_text(out/cfg["outputs"]["cross_cohort_summary"],txt); write_text(out/cfg["outputs"]["concordance_pass"],"status=PASS\nstage=MDV7_CROSS_COHORT_CONCORDANCE\nreplication_language_allowed=NO\n")
    print(txt)
if __name__=="__main__": main()
