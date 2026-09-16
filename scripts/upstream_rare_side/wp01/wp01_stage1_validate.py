#!/usr/bin/env python3
from pathlib import Path
import sys,json,pandas as pd,numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from wp01_common import *
from wp01_lineage import verify_stable

def main():
    root=Path(sys.argv[1]).resolve()
    ma=root/"12_molecular_autism_revision"
    out=ma/"01_anchor_ingest/WP01_raw26_primary_v1"
    verify_stable(ma,out,sha256_file(out/"WP01_preanalysis_lock.json"))
    req=[
      "WP01_GLOBAL_CHECKPOINT_PASS.txt","WP01_ALIGNMENT_PASS.txt","WP01_ALIGNMENT_LOCK.json",
      "WP01_model_frame.tsv.gz","WP01_SOLVER_FULL_REPLAY_PASS.txt","WP01_solver_full_replay_summary.tsv",
      "WP01_solver_full_replay.tsv.gz","WP01_legacy_solver_contract.tsv","WP01_solver_environment.tsv",
      "WP01_RAW26_MDV3_PASS.txt","WP01_raw26_mdv3_all.tsv.gz","WP01_raw26_mdv3_hits.tsv",
      "WP01_raw26_candidate_pools_K50.tsv.gz","WP01_raw26_pilot_sets_1k.tsv.gz",
      "WP01_raw26_pilot_balance_summary.tsv"
    ]
    rows=[]
    def add(x,status,details=""): rows.append({"check_id":x,"status":status,"details":details})
    for f in req:
        p=out/f
        add("OUTPUT::"+f,"PASS" if p.exists() and p.stat().st_size>0 else "FAIL",str(p))

    try:
        sr=pd.read_csv(out/"WP01_solver_full_replay_summary.tsv",sep="\t").iloc[0]
        ok=(int(sr.replay_n)==6671 and int(sr.replay_pass_n)==6671 and int(sr.replay_fail_n)==0
            and int(sr.solver_attempt_mismatch_n)==0 and int(sr.solver_type_mismatch_n)==0
            and int(sr.enriched_mismatch_n)==0 and int(sr.depleted_mismatch_n)==0
            and int(sr.noncertified_n)==0)
        add("SOLVER::full_replay","PASS" if ok else "FAIL",sr.to_json())
    except Exception as e:add("SOLVER::full_replay","FAIL",repr(e))

    try:
        c=pd.read_csv(out/"WP01_legacy_solver_contract.tsv",sep="\t")
        a=c[c.id=="A_ASMEAN_RA05"]
        add("SOLVER::contract","PASS" if len(a)==1 and str(a.iloc[0]["type"])=="AS_mean" else "FAIL",c.to_json())
    except Exception as e:add("SOLVER::contract","FAIL",repr(e))

    try:
        d=pd.read_csv(out/"WP01_raw26_mdv3_all.tsv.gz",sep="\t")
        src=d.groupby("source").size().to_dict()
        fam=(len(d)==6671 and sorted(src.values())==[1221,5450])
        add("MDV3::pathway_family","PASS" if fam else "FAIL",str(src))
        fail=(d.fit_status.astype(str)!="PASS").sum()
        add("MDV3::fit_failures","PASS" if fail==0 else "FAIL",f"fail={fail}")
        add("MDV3::solver_metadata","PASS" if d.solver_attempt.notna().all() and d.solver_type.notna().all() else "FAIL",
            f"attempts={d.solver_attempt.value_counts().to_dict()}; types={d.solver_type.value_counts().to_dict()}")
    except Exception as e:
        add("MDV3::pathway_family","FAIL",repr(e)); add("MDV3::fit_failures","FAIL",repr(e)); add("MDV3::solver_metadata","FAIL",repr(e))

    try:
        pools=pd.read_csv(out/"WP01_raw26_candidate_pools_K50.tsv.gz",sep="\t")
        sizes=pools.groupby("anchor_gene").size()
        add("K50::pool_support","PASS" if len(sizes)==26 and (sizes==50).all() else "FAIL",
            f"anchors={len(sizes)} min={sizes.min()} max={sizes.max()}")
        s=pd.read_csv(out/"WP01_raw26_pilot_sets_1k.tsv.gz",sep="\t")
        add("K50::pilot_unique","PASS" if len(s)==1000 and s.set_sha256.nunique()==1000 else "FAIL",
            f"n={len(s)} unique={s.set_sha256.nunique()}")
        b=pd.read_csv(out/"WP01_raw26_pilot_balance_summary.tsv",sep="\t")
        add("K50::pilot_balance","PASS" if len(b)==3 and (b.status=="PASS").all() else "FAIL",b.to_json())
    except Exception as e:
        add("K50::pool_support","FAIL",repr(e)); add("K50::pilot_unique","FAIL",repr(e)); add("K50::pilot_balance","FAIL",repr(e))

    df=pd.DataFrame(rows)
    df.to_csv(out/"WP01_stage1_checks.tsv",sep="\t",index=False)
    bad=df[df.status!="PASS"]
    for f in ["WP01_STAGE1_HOLD.txt","WP01_STAGE1_PASS.txt","WP01_STOPB_READY.txt","WP01_stage1_lock.json"]:
        try:(out/f).unlink()
        except FileNotFoundError:pass
    if len(bad):
        (out/"WP01_STAGE1_HOLD.txt").write_text("HOLD\n",encoding="utf-8")
        print(df.to_string(index=False)); return 3

    files=req+["WP01_stage1_checks.tsv"]
    obj={
      "stage":"WP01","gate":"B","version":"solver-locked-v1.3.0",
      "preanalysis_lock_sha256":sha256_file(out/"WP01_preanalysis_lock.json"),
      "alignment_lock_sha256":sha256_file(out/"WP01_ALIGNMENT_LOCK.json"),
      "solver_replay_pass_sha256":sha256_file(out/"WP01_SOLVER_FULL_REPLAY_PASS.txt"),
      "raw26_mdv3_pass_sha256":sha256_file(out/"WP01_RAW26_MDV3_PASS.txt"),
      "files":{f:sha256_file(out/f) for f in files},"status":"PASS_CANDIDATE"
    }
    write_json(out/"WP01_stage1_lock.json",obj)
    sh=sha256_file(out/"WP01_stage1_lock.json")
    (out/"WP01_STOPB_READY.txt").write_text(
      f"status=READY_FOR_MANUAL_REVIEW\nstage=WP01\ngate=B\n"
      f"preanalysis_lock_sha256={obj['preanalysis_lock_sha256']}\n"
      f"alignment_lock_sha256={obj['alignment_lock_sha256']}\n"
      f"solver_replay_pass_sha256={obj['solver_replay_pass_sha256']}\n"
      f"raw26_mdv3_pass_sha256={obj['raw26_mdv3_pass_sha256']}\n"
      f"stage1_lock_sha256={sh}\n",encoding="utf-8")
    (out/"WP01_STAGE1_PASS.txt").write_text("PASS\n",encoding="utf-8")
    print((out/"WP01_STOPB_READY.txt").read_text())
    return 0

if __name__=="__main__":
    raise SystemExit(main())
