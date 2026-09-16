#!/usr/bin/env python3
from pathlib import Path
import sys,json,pandas as pd,numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from wp01_common import *
from wp01_lineage import verify_stable

def main():
    root=Path(sys.argv[1]).resolve(); ma=root/"12_molecular_autism_revision"
    out=ma/"01_anchor_ingest/WP01_raw26_primary_v1"
    verify_stable(ma,out,sha256_file(out/"WP01_preanalysis_lock.json"))
    inputs=json.loads((out/"WP01_resolved_inputs.json").read_text())
    rows=[]
    def add(k,s,d=""): rows.append({"check_id":k,"status":s,"details":d})

    for role,info in inputs.items():
        p=Path(info["path"])
        if not p.exists(): add("LEGACY_UNCHANGED::"+role,"FAIL","missing")
        else:
            sh=sha256_file(p)
            add("LEGACY_UNCHANGED::"+role,"PASS" if sh==info["sha256"] else "FAIL",f"sha256={sh}")

    required=[
      "WP01_GLOBAL_CHECKPOINT_PASS.txt","WP01_ALIGNMENT_PASS.txt","WP01_ALIGNMENT_LOCK.json","WP01_model_frame.tsv.gz",
      "WP01_SOLVER_FULL_REPLAY_PASS.txt","WP01_solver_full_replay_summary.tsv",
      "WP01_RAW26_MDV3_PASS.txt","WP01_raw26_mdv3_all.tsv.gz",
      "WP01_raw26_full_sets_10k.tsv.gz","WP01_raw26_K20_sets_10k.tsv.gz",
      "WP01_raw26_full_balance_summary.tsv","WP01_raw26_K20_balance_summary.tsv",
      "WP01_raw26_tier.tsv.gz","WP01_raw26_tier1.tsv","WP01_raw26_vs_core25.tsv.gz","WP01_raw26_vs_core25_summary.tsv"
    ]
    for f in required:
        p=out/f; add("OUTPUT::"+f,"PASS" if p.exists() and p.stat().st_size>0 else "FAIL",str(p))

    try:
        al=json.loads((out/"WP01_ALIGNMENT_LOCK.json").read_text())
        add("ALIGNMENT::model_frame_hash",
            "PASS" if sha256_file(out/"WP01_model_frame.tsv.gz")==al["model_frame_sha256"] else "FAIL",
            sha256_file(out/"WP01_model_frame.tsv.gz"))
    except Exception as e:add("ALIGNMENT::model_frame_hash","FAIL",repr(e))

    try:
        sr=pd.read_csv(out/"WP01_solver_full_replay_summary.tsv",sep="\t").iloc[0]
        ok=(int(sr.replay_n)==6671 and int(sr.replay_pass_n)==6671 and int(sr.replay_fail_n)==0
            and int(sr.solver_attempt_mismatch_n)==0 and int(sr.solver_type_mismatch_n)==0
            and int(sr.enriched_mismatch_n)==0 and int(sr.depleted_mismatch_n)==0
            and int(sr.noncertified_n)==0)
        add("SOLVER::full_replay","PASS" if ok else "FAIL",sr.to_json())
    except Exception as e:add("SOLVER::full_replay","FAIL",repr(e))

    try:
        d=pd.read_csv(out/"WP01_raw26_mdv3_all.tsv.gz",sep="\t")
        src=d.groupby("source").size().to_dict()
        add("MDV3::family","PASS" if len(d)==6671 and sorted(src.values())==[1221,5450] else "FAIL",str(src))
        fail=(d.fit_status.astype(str)!="PASS").sum()
        add("MDV3::fit_fail","PASS" if fail==0 else "FAIL",f"fail={fail}")
    except Exception as e:add("MDV3::family","FAIL",repr(e)); add("MDV3::fit_fail","FAIL",repr(e))

    try:
        pilot=pd.read_csv(out/"WP01_raw26_pilot_sets_1k.tsv.gz",sep="\t")
        full=pd.read_csv(out/"WP01_raw26_full_sets_10k.tsv.gz",sep="\t")
        k20=pd.read_csv(out/"WP01_raw26_K20_sets_10k.tsv.gz",sep="\t")
        add("MDV4::K50_full_unique","PASS" if len(full)==10000 and full.set_sha256.nunique()==10000 else "FAIL",
            f"n={len(full)} unique={full.set_sha256.nunique()}")
        add("MDV4::K20_full_unique","PASS" if len(k20)==10000 and k20.set_sha256.nunique()==10000 else "FAIL",
            f"n={len(k20)} unique={k20.set_sha256.nunique()}")
        same=np.array_equal(pilot.sort_values("replicate_id").set_sha256.to_numpy(),
                            full.sort_values("replicate_id").head(1000).set_sha256.to_numpy())
        add("MDV4::pilot_prefix_reproduced","PASS" if same else "FAIL",str(same))
        b50=pd.read_csv(out/"WP01_raw26_full_balance_summary.tsv",sep="\t")
        b20=pd.read_csv(out/"WP01_raw26_K20_balance_summary.tsv",sep="\t")
        add("MDV4::K50_balance","PASS" if len(b50)==3 and (b50.status=="PASS").all() else "FAIL",b50.to_json())
        add("MDV4::K20_balance","PASS" if len(b20)==3 and (b20.status=="PASS").all() else "FAIL",b20.to_json())
    except Exception as e:
        for k in ["K50_full_unique","K20_full_unique","pilot_prefix_reproduced","K50_balance","K20_balance"]:
            add("MDV4::"+k,"FAIL",repr(e))

    try:
        t=pd.read_csv(out/"WP01_raw26_tier.tsv.gz",sep="\t")
        rule=(t.OR_Firth>1)&(t.q_Firth<0.05)&(t.q_emp<0.05)
        add("TIER::logic","PASS" if np.array_equal(rule.to_numpy(),(t.tier_class=="Tier1").to_numpy()) else "FAIL",
            f"Tier1={int(rule.sum())}")
        src=t.groupby("source").size().to_dict()
        add("TIER::family","PASS" if len(t)==6671 and sorted(src.values())==[1221,5450] else "FAIL",str(src))
    except Exception as e:add("TIER::logic","FAIL",repr(e)); add("TIER::family","FAIL",repr(e))

    df=pd.DataFrame(rows); df.to_csv(out/"WP01_postflight_checks.tsv",sep="\t",index=False)
    bad=df[df.status!="PASS"]
    legacy_changed=int(any((df.check_id.str.startswith("LEGACY_UNCHANGED"))&(df.status!="PASS")))
    status="PASS" if len(bad)==0 else "HOLD"
    summary=f"WP01 validation summary\nstatus={status}\npostflight_failures={len(bad)}\nlegacy_changed={legacy_changed}\n"
    (out/"WP01_validation_summary.txt").write_text(summary,encoding="utf-8")
    lines=[]
    for p in sorted(out.glob("WP01_*")):
        if p.is_file() and p.name!="WP01_checksums.sha256":
            lines.append(f"{sha256_file(p)}  {p.name}")
    (out/"WP01_checksums.sha256").write_text("\n".join(lines)+"\n",encoding="utf-8")
    for m in ["WP01_PASS.txt","WP01_HOLD.txt"]:
        try:(out/m).unlink()
        except FileNotFoundError:pass
    (out/("WP01_PASS.txt" if status=="PASS" else "WP01_HOLD.txt")).write_text(status+"\n",encoding="utf-8")
    print(summary)
    return 0 if status=="PASS" else 3

if __name__=="__main__":
    raise SystemExit(main())
