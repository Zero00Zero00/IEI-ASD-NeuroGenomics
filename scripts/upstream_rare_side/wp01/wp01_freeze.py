#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,sys,pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from wp01_common import *

def main():
    root=Path(sys.argv[1]).resolve()
    ma=root/"12_molecular_autism_revision"; out=ma/"01_anchor_ingest/WP01_raw26_primary_v1"
    if not (out/"WP01_PREFLIGHT_PASS.txt").exists():
        raise RuntimeError("WP01 preflight PASS required")
    contract=read_contract(ma)
    inputs={}
    rels=contract["critical_inputs"]
    for role,rel in rels.items():
        p=root/rel
        if not p.exists(): raise RuntimeError(f"Missing critical input {role}: {p}")
        inputs[role]={"path":str(p),"size":p.stat().st_size,"mtime_ns":p.stat().st_mtime_ns,"sha256":sha256_file(p)}
    # R00 binding
    r00lock=ma/"00_contract/MA_REVISION_PROJECT_LOCK.json"
    inputs["MA_R00_PROJECT_LOCK"]={"path":str(r00lock),"size":r00lock.stat().st_size,"mtime_ns":r00lock.stat().st_mtime_ns,"sha256":sha256_file(r00lock)}
    write_json(out/"WP01_resolved_inputs.json",inputs)

    frame,_=load_structural_frame(root)
    raw=frame[frame.IUIS_flag & frame.SFARI_R0_highconf_flag].copy()
    core=frame[frame.CoreSeed_flag].copy()
    raw[["gene","IUIS_flag","SFARI_R0_highconf_flag","CoreSeed_flag","R0_group"]].to_csv(out/"WP01_raw26_gene_set.tsv",sep="\t",index=False)
    core[["gene","IUIS_flag","SFARI_R0_highconf_flag","CoreSeed_flag","R0_group"]].to_csv(out/"WP01_core25_historical_gene_set.tsv",sep="\t",index=False)
    c4=raw[raw.gene=="C4B"]
    c4[["gene","IUIS_flag","SFARI_R0_highconf_flag","CoreSeed_flag","R0_group"]].to_csv(out/"WP01_C4B_membership.tsv",sep="\t",index=False)

    # code/config hashes
    code_files=[]
    for p in sorted((ma/"workflow/scripts").glob("wp01_*"))+sorted((ma/"bin").glob("*wp01*"))+sorted((ma/"config").glob("wp01_*")):
        if p.is_file():
            code_files.append({"path":str(p),"sha256":sha256_file(p)})
    write_json(out/"WP01_code_lock.json",{"files":code_files})

    resolved_sha=sha256_file(out/"WP01_resolved_inputs.json")
    code_sha=sha256_file(out/"WP01_code_lock.json")
    lock={
      "track":"MOLECULAR_AUTISM_REVISION","stage":"WP01","gate":"A","package_version":PACKAGE_VERSION,
      "ch3_root":str(root),"ma_root":str(ma),"resolved_inputs_sha256":resolved_sha,
      "code_lock_sha256":code_sha,"contract_sha256":sha256_file(ma/"config/wp01_contract.json"),
      "primary_anchor":"raw26","historical_sensitivity":"CoreSeed25",
      "rules":contract["rules"],"seeds":contract["seeds"],
      "data_firewall":"Legacy MDV0-8 and MA-R00 are READ ONLY; WP01 writes only under 12_molecular_autism_revision."
    }
    write_json(out/"WP01_preanalysis_lock.json",lock)
    locksha=sha256_file(out/"WP01_preanalysis_lock.json")
    (out/"WP01_STOPA_READY.txt").write_text(
      "status=READY_FOR_MANUAL_REVIEW\nstage=WP01\ngate=A\n"
      f"preanalysis_lock_sha256={locksha}\nraw26_sha256={sha256_file(out/'WP01_raw26_gene_set.tsv')}\n"
      f"core25_sha256={sha256_file(out/'WP01_core25_historical_gene_set.tsv')}\n",encoding="utf-8")
    print((out/"WP01_STOPA_READY.txt").read_text())
if __name__=="__main__": main()
