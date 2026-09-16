#!/usr/bin/env python3
import argparse, json
from pathlib import Path
import pandas as pd
from mdv7_common import *

def parse_target(df,t):
    for col in ["FULL_NAME","VARIABLE"]:
        if col in df.columns:
            x=df[df[col].astype(str)==t]
            if len(x): return x.iloc[0]
    raise RuntimeError(f"Target {t} absent from MAGMA GSA output")

def family(cfg,out,magma,graw,setfile,scenario,targets,conditional=True):
    sd=out/"gsa"/scenario; sd.mkdir(parents=True,exist_ok=True); pref=sd/"marginal"
    run([magma,"--gene-results",graw,"--set-annot",setfile,"--model",f"direction={cfg['magma']['gene_set_direction']}","--out",pref],str(pref)+".wrapper.log")
    mg=read_magma_table(str(pref)+".gsa.out"); rows=[]
    for t in targets:
        r=parse_target(mg,t); rows.append({"unit_id":t,"NGENES":int(float(r.NGENES)),"BETA":float(r.BETA),"BETA_STD":float(r.get("BETA_STD","nan")),"SE":float(r.SE),"P":float(r.P)})
    marg=pd.DataFrame(rows); marg["q_BH"]=bh(marg.P); cond=pd.DataFrame()
    if conditional and len(targets)>1:
        rows=[]
        for t in targets:
            others=[x for x in targets if x!=t]; tf=sd/f"target_{t}.txt"; tf.write_text(t+"\n"); pref=sd/f"conditional_{t}"; model=[f"analyse=file,{tf}",f"condition-hide={','.join(others)}",f"direction={cfg['magma']['gene_set_direction']}"]
            run([magma,"--gene-results",graw,"--set-annot",setfile,"--model",*model,"--out",pref],str(pref)+".wrapper.log"); cg=read_magma_table(str(pref)+".gsa.out"); r=parse_target(cg,t); rows.append({"unit_id":t,"NGENES":int(float(r.NGENES)),"BETA":float(r.BETA),"BETA_STD":float(r.get("BETA_STD","nan")),"SE":float(r.SE),"P":float(r.P)})
        cond=pd.DataFrame(rows); cond["q_BH"]=bh(cond.P)
    marg.to_csv(sd/"marginal.tsv",sep="\t",index=False)
    if len(cond): cond.to_csv(sd/"conditional.tsv",sep="\t",index=False)
    return marg,cond

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); ap.add_argument("--scenario",required=True,choices=["expanded","compact","w0","w50","mhc","remove_coreseed","alltier1","nested"]); a=ap.parse_args(); cfg=load_cfg(a.config); out=out_dir(cfg); res=json.load(open(out/cfg["outputs"]["resolved_inputs"])); magma=res["magma_bin"]; units=cfg["expected"]["unit_ids"]
    gscenario={"w0":"w0","w50":"w50","mhc":"mhc"}.get(a.scenario,"w10"); graw=out/"magma"/f"SPARK_EUR_{gscenario}.genes.raw"; require_file(graw)
    if a.scenario in ["expanded","w0","w50","mhc"]: setfile=out/cfg["outputs"]["primary_sets"]; targets=units; conditional=True
    elif a.scenario=="compact": setfile=out/cfg["outputs"]["compact_sets"]; targets=units; conditional=True
    elif a.scenario=="remove_coreseed": setfile=out/cfg["outputs"]["remove_coreseed_sets"]; targets=units; conditional=True
    elif a.scenario=="alltier1": setfile=out/cfg["outputs"]["all_tier1_sets"]; targets=units; conditional=False
    else: setfile=out/cfg["outputs"]["nested_sets"]; targets=["D01","M01","M02","M03","M04_EXCLUSIVE","M05","M06"]; conditional=True
    marg,cond=family(cfg,out,magma,graw,setfile,a.scenario,targets,conditional)
    if a.scenario=="expanded":
        m=marg.rename(columns={"NGENES":"NGENES_marginal","BETA":"beta_marginal","BETA_STD":"BETA_STD_marginal","SE":"SE_marginal","P":"P_marginal","q_BH":"q_marginal"}); c=cond.rename(columns={"NGENES":"NGENES_conditional","BETA":"beta_conditional","BETA_STD":"BETA_STD_conditional","SE":"SE_conditional","P":"P_conditional","q_BH":"q_conditional"}); z=m.merge(c,on="unit_id")
        z["primary_support_flag"]=(z.beta_marginal>0)&(z.q_marginal<cfg["primary_inference"]["marginal_family_bh_alpha"]); z["conditional_support_flag"]=(z.beta_conditional>0)&(z.q_conditional<cfg["primary_inference"]["conditional_family_bh_alpha"]); z["interpretation"]="not_detected"; z.loc[z.primary_support_flag,"interpretation"]="independent_SPARK_common_variant_support_with_cross_cohort_heterogeneity"; z["conditional_note"]=""; z.loc[z.unit_id=="M03","conditional_note"]=cfg["interpretation_lock"]["m03_joint_note"]
        m.to_csv(out/cfg["outputs"]["expanded_marginal"],sep="\t",index=False); c.to_csv(out/cfg["outputs"]["expanded_conditional"],sep="\t",index=False); z.to_csv(out/cfg["outputs"]["expanded_summary"],sep="\t",index=False)
        # CoreSeed competitive on same w10 gene results
        csfile=out/cfg["outputs"]["coreseed_set"]; cp=out/"gsa"/"expanded"/"coreseed"; run([magma,"--gene-results",graw,"--set-annot",csfile,"--model",f"direction={cfg['magma']['gene_set_direction']}","--out",cp],str(cp)+".wrapper.log"); cg=read_magma_table(str(cp)+".gsa.out"); r=parse_target(cg,"CoreSeed"); pd.DataFrame([{"set":"CoreSeed","NGENES":int(float(r.NGENES)),"BETA":float(r.BETA),"BETA_STD":float(r.get("BETA_STD","nan")),"SE":float(r.SE),"P":float(r.P)}]).to_csv(out/cfg["outputs"]["coreseed_competitive"],sep="\t",index=False)
        write_text(out/cfg["outputs"]["expanded_pass"],f"status=PASS\nstage=MDV7_PRIMARY_EXPANDED_GSA\nrole=INDEPENDENT_SPARK_FOLLOWUP_NOT_REPLICATION\nunit_n=7\nprimary_support_n={int(z.primary_support_flag.sum())}\nconditionally_supported_n={int(z.conditional_support_flag.sum())}\n")
    elif a.scenario=="compact":
        m=marg.rename(columns={"NGENES":"NGENES_marginal","BETA":"beta_marginal","BETA_STD":"BETA_STD_marginal","SE":"SE_marginal","P":"P_marginal","q_BH":"q_marginal"}); c=cond.rename(columns={"NGENES":"NGENES_conditional","BETA":"beta_conditional","BETA_STD":"BETA_STD_conditional","SE":"SE_conditional","P":"P_conditional","q_BH":"q_conditional"}); z=m.merge(c,on="unit_id"); z["secondary_FDR_flag"]=(z.beta_marginal>0)&(z.q_marginal<cfg["secondary_inference"]["compact_family_bh_alpha"]); z["M02_preidentified_candidate"]=(z.unit_id=="M02"); m.to_csv(out/cfg["outputs"]["compact_marginal"],sep="\t",index=False); c.to_csv(out/cfg["outputs"]["compact_conditional"],sep="\t",index=False); z.to_csv(out/cfg["outputs"]["compact_summary"],sep="\t",index=False); m02=z[z.unit_id=="M02"].iloc[0]; write_text(out/cfg["outputs"]["compact_pass"],f"status=PASS\nstage=MDV7_COMPACT_SECONDARY\nfamily_n=7\nM02_beta={m02.beta_marginal:.12g}\nM02_P={m02.P_marginal:.12g}\nM02_q={m02.q_marginal:.12g}\nM02_secondary_FDR_flag={bool(m02.secondary_FDR_flag)}\nrole=SECONDARY_NEVER_RESCUE_EXPANDED_PRIMARY\n")
    print("MDV7_GSA_DONE",a.scenario)
if __name__=="__main__": main()
