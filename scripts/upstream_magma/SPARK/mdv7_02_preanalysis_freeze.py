#!/usr/bin/env python3
import argparse, json, shutil
from pathlib import Path
import pandas as pd
from mdv7_common import *

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); a=ap.parse_args(); cfg=load_cfg(a.config); root=Path(cfg["project_root"]); out=out_dir(cfg); prov=prov_dir(cfg); require_file(out/cfg["outputs"]["preflight_pass"])
    mapping=[
      ("mdv6_gene_id_map","gene_id_map"),("mdv6_gene_locations","gene_locations"),("mdv6_primary_sets","primary_sets"),("mdv6_compact_sets","compact_sets"),("mdv6_remove_coreseed_sets","remove_coreseed_sets"),("mdv6_all_tier1_sets","all_tier1_sets"),("mdv6_nested_sets","nested_sets"),("mdv6_coreseed_set","coreseed_set")]
    rows=[]
    for src_key,out_key in mapping:
        src=root_path(cfg,cfg["upstream"][src_key]); dst=out/cfg["outputs"][out_key]; require_file(src); shutil.copy2(src,dst); same=sha256_file(src)==sha256_file(dst); rows.append({"source_key":src_key,"source_path":str(src.relative_to(root)),"mdv7_path":str(dst.relative_to(root)),"source_sha256":sha256_file(src),"mdv7_sha256":sha256_file(dst),"byte_identical":same})
        if not same: raise RuntimeError(f"Byte-identical reuse failed: {src_key}")
    pd.DataFrame(rows).to_csv(out/cfg["outputs"]["freeze_manifest"],sep="\t",index=False)
    # verify expected set counts from copied files
    def read_sets(path):
        d={}
        for line in Path(path).read_text().splitlines():
            a=line.split(); d[a[0]]=set(a[1:])
        return d
    exp=read_sets(out/cfg["outputs"]["primary_sets"]); comp=read_sets(out/cfg["outputs"]["compact_sets"])
    if sorted(exp)!=cfg["expected"]["unit_ids"] or sorted(comp)!=cfg["expected"]["unit_ids"]: raise RuntimeError("Copied set IDs differ from frozen 7-unit family")
    for u,n in cfg["expected"]["expanded_counts"].items():
        if len(exp[u])!=n: raise RuntimeError(f"{u} Expanded count drift: {len(exp[u])} != {n}")
    for u,n in cfg["expected"]["compact_counts"].items():
        if len(comp[u])!=n: raise RuntimeError(f"{u} Compact count drift: {len(comp[u])} != {n}")
    # code/config/interpretation freeze
    pkg=Path(a.config).parent.parent; code=[Path(a.config),pkg/"Snakefile.mdv7",pkg/"MDV7_PREANALYSIS_DECISION_LOCK.md",pkg/"workflow/envs/mdv7.yaml",*sorted((pkg/"workflow/scripts").glob("mdv7_*.py"))]
    codehash={str(p.relative_to(pkg)):sha256_file(p) for p in code if p.exists()}
    scientific={r["mdv7_path"]:r["mdv7_sha256"] for r in rows}
    lock={"version":cfg["version"],"stage":cfg["stage"],"association_rows_read":"NO","SPARK_role":"INDEPENDENT_COMMON_VARIANT_FOLLOWUP_NOT_REPLICATION","SPARK_primary":"SPARK_ONLY_EUR","primary_gene_window_kb":10,"primary_gene_sets":"Domain-Expanded","primary_units":cfg["expected"]["unit_ids"],"secondary_compact_family":cfg["expected"]["unit_ids"],"M02_compact_preidentified_secondary_candidate":True,"forbidden_meta_files":[cfg["spark"]["excluded_meta_eur"],cfg["spark"]["excluded_meta_all"]],"multiancestry_inferential_use":False,"harmonization_strategy":cfg["harmonization"]["strategy"],"interpretation_lock":cfg["interpretation_lock"],"scientific_files":scientific,"code_and_config_sha256":codehash,"input_lock_sha256":sha256_file(out/cfg["outputs"]["input_lock"])}
    write_json(out/cfg["outputs"]["preanalysis_lock"],lock)
    # copy decision lock into provenance before association rows are read
    shutil.copy2(pkg/"MDV7_PREANALYSIS_DECISION_LOCK.md",prov/"MDV7_PREANALYSIS_DECISION_LOCK.md")
    write_text(out/cfg["outputs"]["preanalysis_pass"],"status=PASS\ntechnical_status=PASS\nstage=MDV7_PREANALYSIS_FREEZE\nassociation_rows_read=NO\nprimary_dataset=SPARK_ONLY_EUR\nprimary_window_kb=10\nprimary_set=Domain-Expanded\nsecondary_compact_family=7_units\nM02_compact_candidate=PREIDENTIFIED_SECONDARY\nforbidden_PGC_overlap_meta_files=YES\n")
    print("MDV7_PREANALYSIS_FREEZE_PASS")
if __name__=="__main__": main()
