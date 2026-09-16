#!/usr/bin/env python3
import argparse, gzip, json, os, sys
from pathlib import Path
import pandas as pd
from mdv6_common import *

def add(rows,category,check,status,observed,expected="",detail=""):
    rows.append(dict(category=category,check=check,status=status,observed=observed,expected=expected,detail=detail))
    if status=="FAIL": print(f"[FAIL] {category}::{check}: {detail or observed}")

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); a=ap.parse_args()
    cfg=load_cfg(a.config); root=Path(cfg["project_root"]); out=out_dir(cfg); prov=prov_dir(cfg); rows=[]
    # required upstream PASS and checksum manifests
    for st in ["mdv0","mdv1","mdv2","mdv3","mdv4","mdv5"]:
        pp=root_path(cfg,cfg["upstream"][f"{st}_pass"])
        add(rows,"upstream",f"{st}_pass","PASS" if pp.exists() and pp.stat().st_size else "FAIL",str(pp),"existing non-empty")
        cm=root_path(cfg,cfg["upstream"][f"{st}_checksums"])
        vr=verify_sha256_manifest(cm,root)
        fail=sum(x["status"]!="PASS" for x in vr)
        add(rows,"checksums",f"{st}_checksums","PASS" if fail==0 and len(vr)>0 else "FAIL",f"verified={len(vr)-fail}; failed={fail}","0 failed")
    # mandatory files
    keys=["mdv4_pseudo_sets","mdv5_analysis_lock","mdv5_domain_membership","gene_universe","mdv5_biological_labels","mdv5_final_biological_audit"]
    for k in keys:
        p=root_path(cfg,cfg["upstream"][k]); add(rows,"input",k,"PASS" if p.exists() and p.stat().st_size else "FAIL",str(p),"existing non-empty")
    # universe
    up=root_path(cfg,cfg["upstream"]["gene_universe"])
    if up.exists():
        u=pd.read_csv(up,sep="\t")
        add(rows,"universe","n","PASS" if len(u)==cfg["expected"]["universe_n"] else "FAIL",len(u),cfg["expected"]["universe_n"])
        cs=int(pd.to_numeric(u["CoreSeed_flag"],errors="coerce").fillna(0).sum())
        add(rows,"universe","coreseed_n","PASS" if cs==cfg["expected"]["coreseed_n"] else "FAIL",cs,cfg["expected"]["coreseed_n"])
        add(rows,"universe","unique_hgnc","PASS" if u["HGNC_id"].nunique()==len(u) else "FAIL",u["HGNC_id"].nunique(),len(u))
    # MDV5 units
    mp=root_path(cfg,cfg["upstream"]["mdv5_domain_membership"])
    if mp.exists():
        d=pd.read_csv(mp,sep="\t")
        ids=sorted(d["domain_id"].dropna().unique().tolist()); exp=cfg["expected"]["unit_ids"]
        add(rows,"mdv5","unit_ids","PASS" if ids==exp else "FAIL",";".join(ids),";".join(exp))
        for uid,n in cfg["expected"]["expanded_counts"].items():
            got=int(d.loc[(d.domain_id==uid)&(d.Expanded_flag==1),"HGNC_id"].nunique())
            add(rows,"mdv5",f"{uid}_expanded","PASS" if got==n else "FAIL",got,n)
        for uid,n in cfg["expected"]["compact_counts"].items():
            got=int(d.loc[(d.domain_id==uid)&(d.Compact_flag==1),"HGNC_id"].nunique())
            add(rows,"mdv5",f"{uid}_compact","PASS" if got==n else "FAIL",got,n)
    # labels exact
    lp=root_path(cfg,cfg["upstream"]["mdv5_biological_labels"])
    if lp.exists():
        lab=pd.read_csv(lp,sep="\t")
        ids=lab["unit_id"].tolist()
        add(rows,"mdv5","biological_label_freeze","PASS" if ids==cfg["expected"]["unit_ids"] else "FAIL",";".join(ids),";".join(cfg["expected"]["unit_ids"]))
    # PGC provider file identity/header (does not inspect association distribution)
    pgc=root_path(cfg,cfg["pgc2019"]["gwas"]); readme=root_path(cfg,cfg["pgc2019"]["readme"])
    try:
        require_file(pgc,"PGC2019 GWAS"); require_file(readme,"PGC README")
        with gzip.open(pgc,"rt",encoding="utf-8",errors="replace") as f: hdr=f.readline().strip().split()
        expected=cfg["pgc2019"]["expected_header"]
        add(rows,"pgc","header","PASS" if all(x in hdr for x in expected) else "FAIL",";".join(hdr),";".join(expected))
        with open(readme,"rb") as f: pdf=f.read(5)==b"%PDF-"
        add(rows,"pgc","readme_pdf","PASS" if pdf else "FAIL",pdf,True)
        pgc_sha=sha256_file(pgc)
    except Exception as e:
        add(rows,"pgc","identity","FAIL","","",str(e)); pgc_sha=""
    # 1KG
    try:
        kdir=root_path(cfg,cfg["ld_reference"]["extracted_dir"])
        prefix=cfg["ld_reference"].get("prefix","AUTO")
        if prefix=="AUTO": prefix=discover_1kg_prefix(kdir)
        elif not Path(prefix).is_absolute(): prefix=str(root/Path(prefix))
        for ext in ["bed","bim","fam"]: require_file(prefix+"."+ext,f"1KG {ext}")
        add(rows,"ld","prefix","PASS",prefix,"BED/BIM/FAM")
    except Exception as e:
        prefix=""; add(rows,"ld","prefix","FAIL","","",str(e))
    # MAGMA
    try:
        magma=resolve_magma(cfg["magma"].get("binary","AUTO")); ver=magma_version(magma)
        ok=cfg["magma"]["required_version"] in ver
        add(rows,"software","magma_version","PASS" if ok else "FAIL",ver,cfg["magma"]["required_version"])
        magma_sha=sha256_file(magma)
    except Exception as e:
        magma=""; ver=""; magma_sha=""; add(rows,"software","magma_version","FAIL","","",str(e))
    # No SPARK declared inputs
    declared=json.dumps({"upstream":cfg["upstream"],"pgc2019":cfg["pgc2019"],"ld_reference":cfg["ld_reference"]}).lower()
    bad="v1_spark" in declared or "/spark" in declared or "spark_" in declared
    add(rows,"separation","spark_not_declared","PASS" if not bad else "FAIL",not bad,True)
    # report and lock
    report=out/cfg["outputs"]["preflight_report"]; pd.DataFrame(rows).to_csv(report,sep="\t",index=False)
    fails=sum(r["status"]=="FAIL" for r in rows)
    resolved={"project_root":str(root),"pgc_gwas":str(pgc),"pgc_readme":str(readme),"pgc_sha256":pgc_sha,"ld_prefix":prefix,"magma_bin":magma,"magma_version":ver,"magma_sha256":magma_sha,"association_results_read":"NO"}
    write_json(out/cfg["outputs"]["resolved_inputs"],resolved)
    lock={"version":cfg["version"],"stage":cfg["stage"],"association_results_read":"NO","config_sha256":sha256_file(a.config),"inputs":{}}
    for k in keys:
        p=root_path(cfg,cfg["upstream"][k]); lock["inputs"][k]={"path":str(p),"sha256":sha256_file(p) if p.exists() and p.is_file() else None}
    lock["inputs"]["pgc2019"]={"path":str(pgc),"sha256":pgc_sha,"build":cfg["pgc2019"]["build"]}
    lock["inputs"]["ld_prefix"]={"prefix":prefix,"bim_sha256":sha256_file(prefix+".bim") if prefix else None}
    lock["software"]={"magma_bin":magma,"version":ver,"sha256":magma_sha}
    write_json(out/cfg["outputs"]["input_lock"],lock)
    if fails: raise SystemExit(f"MDV6 preflight FAIL: {fails} checks failed; see {report}")
    write_text(out/cfg["outputs"]["preflight_pass"],f"status=PASS\ntechnical_status=PASS\nstage=MDV6_PREFLIGHT\nassociation_results_read=NO\npgc_build={cfg['pgc2019']['build']}\nmagma_version={ver}\n")
    print("MDV6_PREFLIGHT_PASS")
if __name__=="__main__": main()
