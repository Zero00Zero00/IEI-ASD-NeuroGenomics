#!/usr/bin/env python3
import argparse, json, os, tempfile
from pathlib import Path
import pandas as pd
from mdv6_common import *

def parse_target(df,target):
    # prefer FULL_NAME, fallback VARIABLE; strip MAGMA truncation only as fallback
    if "FULL_NAME" in df.columns:
        x=df[df.FULL_NAME.astype(str)==target]
        if len(x): return x.iloc[0]
    if "VARIABLE" in df.columns:
        x=df[df.VARIABLE.astype(str)==target]
        if len(x): return x.iloc[0]
    raise RuntimeError(f"Target {target} not found in MAGMA GSA output")

def run_family(cfg,out,magma,graw,setfile,scenario,targets,conditional=True):
    sdir=out/"gsa"/scenario; sdir.mkdir(parents=True,exist_ok=True)
    pref=sdir/"marginal"
    run([magma,"--gene-results",graw,"--set-annot",setfile,"--model",f"direction={cfg['magma']['gene_set_direction']}","--out",pref],str(pref)+".wrapper.log")
    mg=read_magma_table(str(pref)+".gsa.out")
    mrows=[]
    for t in targets:
        r=parse_target(mg,t)
        mrows.append({"unit_id":t,"NGENES":int(float(r["NGENES"])),"BETA":float(r["BETA"]),"BETA_STD":float(r.get("BETA_STD",'nan')),"SE":float(r["SE"]),"P":float(r["P"])})
    marg=pd.DataFrame(mrows); marg["q_BH"]=bh(marg.P)
    cond=pd.DataFrame()
    if conditional and len(targets)>1:
        crows=[]
        for t in targets:
            others=[x for x in targets if x!=t]
            tf=sdir/f"target_{t}.txt"; tf.write_text(t+"\n",encoding="utf-8")
            cp=sdir/f"conditional_{t}"
            model=[f"analyse=file,{tf}",f"condition-hide={','.join(others)}",f"direction={cfg['magma']['gene_set_direction']}"]
            run([magma,"--gene-results",graw,"--set-annot",setfile,"--model",*model,"--out",cp],str(cp)+".wrapper.log")
            cg=read_magma_table(str(cp)+".gsa.out"); r=parse_target(cg,t)
            crows.append({"unit_id":t,"NGENES":int(float(r["NGENES"])),"BETA":float(r["BETA"]),"BETA_STD":float(r.get("BETA_STD",'nan')),"SE":float(r["SE"]),"P":float(r["P"])})
        cond=pd.DataFrame(crows); cond["q_BH"]=bh(cond.P)
    return marg,cond

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); ap.add_argument("--scenario",required=True,choices=["primary","w0","w50","mhc","compact","remove_coreseed","alltier1","nested"]); a=ap.parse_args()
    cfg=load_cfg(a.config); out=out_dir(cfg); resolved=json.load(open(out/cfg["outputs"]["resolved_inputs"])); magma=resolved["magma_bin"]
    units=cfg["expected"]["unit_ids"]
    if a.scenario in ["primary","compact","remove_coreseed","alltier1","nested"]: gscenario="w10"
    else: gscenario=a.scenario
    graw=out/"magma"/f"PGC2019_{gscenario}.genes.raw"; require_file(graw)
    if a.scenario in ["primary","w0","w50","mhc"]: setfile=out/cfg["outputs"]["primary_sets"]; targets=units; conditional=True
    elif a.scenario=="compact": setfile=out/cfg["outputs"]["compact_sets"]; targets=units; conditional=True
    elif a.scenario=="remove_coreseed": setfile=out/cfg["outputs"]["remove_coreseed_sets"]; targets=units; conditional=True
    elif a.scenario=="alltier1": setfile=out/cfg["outputs"]["all_tier1_sets"]; targets=units; conditional=False
    else:
        setfile=out/cfg["outputs"]["nested_sets"]; targets=["D01","M01","M02","M03","M04_EXCLUSIVE","M05","M06"]; conditional=True
    marg,cond=run_family(cfg,out,magma,graw,setfile,a.scenario,targets,conditional)
    sdir=out/"gsa"/a.scenario
    marg.to_csv(sdir/"marginal.tsv",sep="\t",index=False)
    if len(cond): cond.to_csv(sdir/"conditional.tsv",sep="\t",index=False)
    # CoreSeed competitive only for primary/window/MHC scenarios using combined one-set file
    if a.scenario in ["primary","w0","w50","mhc"]:
        csfile=out/cfg["outputs"]["coreseed_set"]; cp=sdir/"coreseed"
        run([magma,"--gene-results",graw,"--set-annot",csfile,"--model",f"direction={cfg['magma']['gene_set_direction']}","--out",cp],str(cp)+".wrapper.log")
        cs=read_magma_table(str(cp)+".gsa.out"); r=parse_target(cs,"CoreSeed")
        pd.DataFrame([{"set":"CoreSeed","NGENES":int(float(r["NGENES"])),"BETA":float(r["BETA"]),"BETA_STD":float(r.get("BETA_STD",'nan')),"SE":float(r["SE"]),"P":float(r["P"])}]).to_csv(sdir/"coreseed.tsv",sep="\t",index=False)
    if a.scenario=="primary":
        # Stable public output names + interpretive classes.
        marg2=marg.rename(columns={"BETA":"beta_marginal","SE":"SE_marginal","P":"P_marginal","q_BH":"q_marginal","NGENES":"NGENES_marginal","BETA_STD":"BETA_STD_marginal"})
        cond2=cond.rename(columns={"BETA":"beta_conditional","SE":"SE_conditional","P":"P_conditional","q_BH":"q_conditional","NGENES":"NGENES_conditional","BETA_STD":"BETA_STD_conditional"})
        z=marg2.merge(cond2,on="unit_id",how="left")
        z["primary_transfer_flag"]=(z.beta_marginal>0)&(z.q_marginal<cfg["primary_inference"]["marginal_family_bh_alpha"])
        z["conditional_retention_flag"]=(z.beta_conditional>0)&(z.q_conditional<cfg["primary_inference"]["conditional_family_bh_alpha"])
        z["interpretation"]="not_detected"
        z.loc[z.primary_transfer_flag & ~z.conditional_retention_flag,"interpretation"]="marginal_only_correspondence"
        z.loc[z.primary_transfer_flag & z.conditional_retention_flag,"interpretation"]="conditionally_retained_correspondence"
        # M03 nested warning is always descriptive, never changes p/q.
        z["conditional_note"]=""
        z.loc[z.unit_id=="M03","conditional_note"]=cfg["interpretation_lock"]["m03_joint_note"]
        marg2.to_csv(out/cfg["outputs"]["primary_marginal"],sep="\t",index=False)
        cond2.to_csv(out/cfg["outputs"]["primary_conditional"],sep="\t",index=False)
        cs=pd.read_csv(sdir/"coreseed.tsv",sep="\t"); cs.to_csv(out/cfg["outputs"]["coreseed_competitive"],sep="\t",index=False)
        z.to_csv(out/cfg["outputs"]["primary_gsa_summary"],sep="\t",index=False)
        write_text(out/cfg["outputs"]["primary_gsa_pass"],f"status=PASS\nstage=MDV6_PRIMARY_GSA\nassociation_results_read=YES\nprimary_window_kb=10\nprimary_set=Domain-Expanded\nunit_n=7\nprimary_transfer_n={int(z.primary_transfer_flag.sum())}\nconditionally_retained_n={int(z.conditional_retention_flag.sum())}\npgc_role=PRIMARY_TEST_NOT_REPLICATION\n")
    print("GSA_DONE",a.scenario)
if __name__=="__main__": main()
