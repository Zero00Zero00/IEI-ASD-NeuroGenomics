#!/usr/bin/env python3
import argparse, gzip, json, re
from pathlib import Path
import pandas as pd
from mdv7_common import *

def add(rows,cat,check,status,obs="",exp="",detail=""):
    rows.append(dict(category=cat,check=check,status=status,observed=obs,expected=exp,detail=detail))
    if status=="FAIL": print(f"[FAIL] {cat}::{check}: {detail or obs}")

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); a=ap.parse_args(); cfg=load_cfg(a.config); root=Path(cfg["project_root"]); out=out_dir(cfg); rows=[]
    for st in ["mdv0","mdv1","mdv2","mdv3","mdv4","mdv5","mdv6"]:
        pp=root_path(cfg,cfg["upstream"][f"{st}_pass"]); add(rows,"upstream",f"{st}_pass","PASS" if pp.exists() and pp.stat().st_size else "FAIL",str(pp),"existing non-empty")
        vr=verify_sha256_manifest(root_path(cfg,cfg["upstream"][f"{st}_checksums"]),root); bad=sum(x["status"]!="PASS" for x in vr); add(rows,"checksums",f"{st}_checksums","PASS" if bad==0 and len(vr)>0 else "FAIL",f"verified={len(vr)-bad};failed={bad}","0 failed")
    must=["mdv6_analysis_lock","mdv6_preanalysis_lock","mdv6_gene_id_map","mdv6_gene_locations","mdv6_primary_sets","mdv6_compact_sets","mdv6_remove_coreseed_sets","mdv6_all_tier1_sets","mdv6_nested_sets","mdv6_coreseed_set","mdv6_primary_marginal","mdv6_primary_conditional","mdv6_primary_summary","mdv6_sensitivity_summary","mdv6_coreseed_empirical","mdv4_pseudo_sets","gene_universe","mdv5_biological_labels"]
    for k in must:
        p=root_path(cfg,cfg["upstream"][k]); add(rows,"input",k,"PASS" if p.exists() and p.stat().st_size else "FAIL",str(p),"existing non-empty")
    # verify MDV6 result state that motivates MDV7 interpretation
    try:
        ps=pd.read_csv(root_path(cfg,cfg["upstream"]["mdv6_primary_summary"]),sep="\t"); n=int(pd.to_numeric(ps["primary_transfer_flag"],errors="coerce").fillna(0).astype(bool).sum())
        add(rows,"mdv6_state","pgc_primary_transfer_n","PASS" if n==cfg["expected"]["mdv6_pgc_primary_transfer_n"] else "FAIL",n,cfg["expected"]["mdv6_pgc_primary_transfer_n"])
        ss=pd.read_csv(root_path(cfg,cfg["upstream"]["mdv6_sensitivity_summary"]),sep="\t")
        m=ss[(ss.scenario=="compact")&(ss.model=="marginal")&(ss.unit_id=="M02")]
        ok=len(m)==1 and float(m.iloc[0].BETA)>0 and float(m.iloc[0].q_BH)<0.05
        add(rows,"mdv6_state","M02_compact_secondary_positive","PASS" if ok else "FAIL",ok,True)
    except Exception as e: add(rows,"mdv6_state","result_state","FAIL","","",str(e))
    # SPARK files: identity/header only, no association rows parsed
    spark_keys=["primary_eur","multiancestry_sensitivity_candidate","excluded_meta_eur","excluded_meta_all","readme","sumstats_info","repo_commit","repo_method"]
    for k in spark_keys:
        p=root_path(cfg,cfg["spark"][k]); add(rows,"spark_identity",k,"PASS" if p.exists() and p.stat().st_size else "FAIL",str(p),"existing non-empty")
    sp=root_path(cfg,cfg["spark"]["primary_eur"])
    try:
        with gzip.open(sp,"rt",encoding="utf-8",errors="strict") as f: hdr=f.readline().strip().split()
        exp=cfg["spark"]["expected_header"]; add(rows,"spark_identity","primary_header","PASS" if all(x in hdr for x in exp) else "FAIL",";".join(hdr),";".join(exp))
        spark_sha=sha256_file(sp)
    except Exception as e: spark_sha=""; add(rows,"spark_identity","primary_header","FAIL","","",str(e))
    try:
        commit=root_path(cfg,cfg["spark"]["repo_commit"]).read_text().strip(); add(rows,"spark_identity","repo_commit_40hex","PASS" if re.fullmatch(r"[0-9a-fA-F]{40}",commit) else "FAIL",commit,"40-hex commit")
    except Exception as e: commit=""; add(rows,"spark_identity","repo_commit_40hex","FAIL","","",str(e))
    # hard exclusion of overlapping meta files from declared inferential inputs
    declared=json.dumps({"primary":cfg["spark"]["primary_eur"],"ld":cfg["ld_reference"],"magma":cfg["magma"]}).lower()
    overlap_names=[Path(cfg["spark"]["excluded_meta_eur"]).name.lower(),Path(cfg["spark"]["excluded_meta_all"]).name.lower()]
    add(rows,"independence","overlap_meta_not_declared","PASS" if not any(x in declared for x in overlap_names) else "FAIL",True,True)
    add(rows,"independence","primary_is_exact_SPARK_only_EUR","PASS" if Path(cfg["spark"]["primary_eur"]).name=="SPARK_update_model1QC_EUR_only.tsv.gz" else "FAIL",Path(cfg["spark"]["primary_eur"]).name,"SPARK_update_model1QC_EUR_only.tsv.gz")
    add(rows,"independence","multiancestry_inferential_use_false","PASS" if cfg["multiancestry"]["inferential_use"] is False else "FAIL",cfg["multiancestry"]["inferential_use"],False)
    # liftover resources
    for k in ["binary","chain"]:
        p=root_path(cfg,cfg["liftover"][k]); ok=p.exists() and p.stat().st_size>0 and (k!="binary" or os.access(p,os.X_OK)); add(rows,"liftover",k,"PASS" if ok else "FAIL",str(p),"present"+(" and executable" if k=="binary" else ""))
    # LD + MAGMA
    try:
        kdir=root_path(cfg,cfg["ld_reference"]["extracted_dir"]); prefix=cfg["ld_reference"].get("prefix","AUTO"); prefix=discover_1kg_prefix(kdir) if prefix=="AUTO" else str(root_path(cfg,prefix))
        for ext in ["bed","bim","fam"]: require_file(prefix+"."+ext)
        add(rows,"ld","prefix","PASS",prefix,"BED/BIM/FAM")
    except Exception as e: prefix=""; add(rows,"ld","prefix","FAIL","","",str(e))
    try:
        magma=resolve_magma(cfg["magma"].get("binary","AUTO")); ver=magma_version(magma); add(rows,"software","magma","PASS" if cfg["magma"]["required_version"] in ver else "FAIL",ver,cfg["magma"]["required_version"]); magma_sha=sha256_file(magma)
    except Exception as e: magma=""; ver=""; magma_sha=""; add(rows,"software","magma","FAIL","","",str(e))
    report=out/cfg["outputs"]["preflight_report"]; pd.DataFrame(rows).to_csv(report,sep="\t",index=False); fails=sum(x["status"]=="FAIL" for x in rows)
    resolved={"project_root":str(root),"spark_primary":str(sp),"spark_primary_sha256":spark_sha,"spark_commit":commit,"ld_prefix":prefix,"magma_bin":magma,"magma_version":ver,"magma_sha256":magma_sha,"liftover_bin":str(root_path(cfg,cfg["liftover"]["binary"])),"liftover_chain":str(root_path(cfg,cfg["liftover"]["chain"])),"association_rows_read":"NO"}
    write_json(out/cfg["outputs"]["resolved_inputs"],resolved)
    lock={"version":cfg["version"],"stage":cfg["stage"],"association_rows_read":"NO","config_sha256":sha256_file(a.config),"spark_primary":{"path":str(sp),"sha256":spark_sha,"role":"SPARK_ONLY_EUR","raw_build":cfg["spark"]["raw_build"]},"forbidden_meta_files":[cfg["spark"]["excluded_meta_eur"],cfg["spark"]["excluded_meta_all"]],"ld_prefix":prefix,"magma":{"path":magma,"version":ver,"sha256":magma_sha},"liftover":{"binary":resolved["liftover_bin"],"binary_sha256":sha256_file(resolved["liftover_bin"]) if Path(resolved["liftover_bin"]).exists() else None,"chain":resolved["liftover_chain"],"chain_sha256":sha256_file(resolved["liftover_chain"]) if Path(resolved["liftover_chain"]).exists() else None}}
    write_json(out/cfg["outputs"]["input_lock"],lock)
    if fails: raise SystemExit(f"MDV7 preflight FAIL: {fails} checks failed; see {report}")
    write_text(out/cfg["outputs"]["preflight_pass"],f"status=PASS\ntechnical_status=PASS\nstage=MDV7_PREFLIGHT\nassociation_rows_read=NO\nSPARK_role=INDEPENDENT_FOLLOWUP\nSPARK_primary=SPARK_ONLY_EUR\nSPARK_raw_build=GRCh38\ntarget_build=GRCh37\nmagma_version={ver}\n")
    print("MDV7_PREFLIGHT_PASS")
if __name__=="__main__":
    import os
    main()
