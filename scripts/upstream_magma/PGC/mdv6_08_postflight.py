#!/usr/bin/env python3
import argparse, json, os, platform, sys
from pathlib import Path
import pandas as pd, numpy as np, scipy, yaml
from mdv6_common import *

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); a=ap.parse_args(); cfg=load_cfg(a.config); root=Path(cfg["project_root"]); out=out_dir(cfg)
    checks=[]
    def ck(name,ok,obs="",exp=""):
        checks.append({"check":name,"status":"PASS" if ok else "FAIL","observed":obs,"expected":exp})
    for k in ["preflight_pass","preanalysis_pass","pgc_qc_pass","primary_gene_pass","primary_gsa_pass","coreseed_empirical_pass","sensitivity_pass"]:
        p=out/cfg["outputs"][k]; ck(k,p.exists() and p.stat().st_size>0,str(p),"existing non-empty")
    # Revalidate upstream checksum chain again at terminal boundary.
    for st in ["mdv0","mdv1","mdv2","mdv3","mdv4","mdv5"]:
        vr=verify_sha256_manifest(root_path(cfg,cfg["upstream"][f"{st}_checksums"]),root); bad=sum(x["status"]!="PASS" for x in vr); ck(f"{st}_checksums_terminal",bad==0,f"failed={bad}","0")
    # primary results cardinality + BH recomputation
    p=pd.read_csv(out/cfg["outputs"]["primary_gsa_summary"],sep="\t"); ck("primary_unit_n",len(p)==7,len(p),7)
    ck("primary_unit_ids",sorted(p.unit_id.tolist())==cfg["expected"]["unit_ids"],";".join(sorted(p.unit_id)),";".join(cfg["expected"]["unit_ids"]))
    q1=bh(p.P_marginal); q2=bh(p.P_conditional)
    ck("marginal_BH_recalc",np.allclose(q1,p.q_marginal,rtol=0,atol=1e-12,equal_nan=True))
    ck("conditional_BH_recalc",np.allclose(q2,p.q_conditional,rtol=0,atol=1e-12,equal_nan=True))
    ck("finite_primary",np.isfinite(p[["beta_marginal","SE_marginal","P_marginal","q_marginal","beta_conditional","SE_conditional","P_conditional","q_conditional"]].to_numpy(float)).all())
    # nested freeze still exact
    dm=pd.read_csv(root_path(cfg,cfg["upstream"]["mdv5_domain_membership"]),sep="\t")
    S3=set(dm.loc[(dm.domain_id=="M03")&(dm.Expanded_flag==1),"HGNC_id"]); S4=set(dm.loc[(dm.domain_id=="M04")&(dm.Expanded_flag==1),"HGNC_id"])
    ck("M03_nested_M04_terminal",S3<=S4,len(S3&S4),len(S3))
    # PGC role / no Spark in resolved input
    resolved=json.load(open(out/cfg["outputs"]["resolved_inputs"])); ck("PGC_role_primary_not_replication",cfg["interpretation_lock"]["pgc_role"]=="primary_common_variant_test_not_replication")
    ck("SPARK_not_in_resolved_inputs","spark" not in json.dumps(resolved).lower())
    # CoreSeed null exact fixed source
    ce=pd.read_csv(out/cfg["outputs"]["coreseed_empirical"],sep="\t"); ck("CoreSeed_empirical_two_stats",set(ce.statistic)=={"mean_Z","median_Z"},";".join(ce.statistic))
    # sensitivity scenarios
    ss=pd.read_csv(out/cfg["outputs"]["sensitivity_summary"],sep="\t"); expected=set(["w0","w50","mhc","compact","remove_coreseed","alltier1","nested"]); ck("sensitivity_scenarios",set(ss.scenario)==expected,";".join(sorted(set(ss.scenario))),";".join(sorted(expected)))
    cdf=pd.DataFrame(checks); cdf.to_csv(out/cfg["outputs"]["postflight_checks"],sep="\t",index=False); failed=int((cdf.status=="FAIL").sum())
    # session
    res=json.load(open(out/cfg["outputs"]["resolved_inputs"]))
    write_text(out/cfg["outputs"]["session_versions"],f"python={sys.version}\npandas={pd.__version__}\nnumpy={np.__version__}\nscipy={scipy.__version__}\npyyaml={yaml.__version__}\nplatform={platform.platform()}\nmagma={res.get('magma_version','')}\n")
    if failed: raise RuntimeError(f"MDV6 postflight failed {failed} checks")
    transfer=p[p.primary_transfer_flag==True]; cond=p[p.conditional_retention_flag==True]
    summary=("MDV-6 PGC2019 COMMON-VARIANT CORRESPONDENCE — VALIDATION SUMMARY\n"+"="*70+"\n"+
             f"version={cfg['version']}\nprimary_window_kb=10\nprimary_gene_sets=Domain-Expanded\nprimary_unit_n=7\nprimary_transfer_n={len(transfer)}\nconditionally_retained_n={len(cond)}\nprimary_transfer_units={';'.join(transfer.unit_id) if len(transfer) else 'NONE'}\nconditional_units={';'.join(cond.unit_id) if len(cond) else 'NONE'}\nPGC_role=PRIMARY_COMMON_VARIANT_TEST_NOT_REPLICATION\nSPARK_inputs_used=NO\ntechnical_checks_passed_n={len(cdf)}\ntechnical_checks_failed_n=0\n")
    write_text(out/cfg["outputs"]["validation_summary"],summary)
    lock={"version":cfg["version"],"stage":cfg["stage"],"association_results_read":"YES","PGC_role":"PRIMARY_TEST_NOT_REPLICATION","SPARK_inputs_used":"NO","primary_window_kb":10,"primary_set":"Domain-Expanded","primary_transfer_units":transfer.unit_id.tolist(),"conditionally_retained_units":cond.unit_id.tolist(),"input_lock_sha256":sha256_file(out/cfg["outputs"]["input_lock"]),"preanalysis_lock_sha256":sha256_file(out/cfg["outputs"]["preanalysis_lock"]),"interpretation_lock":cfg["interpretation_lock"]}
    write_json(out/cfg["outputs"]["analysis_lock"],lock)
    # checksums of all canonical outputs except checksums and PASS themselves
    exclude={cfg["outputs"]["checksums"],cfg["outputs"]["pass"]}
    files=[]
    for pth in sorted(out.rglob("*")):
        if pth.is_file() and pth.name not in exclude:
            # Keep MAGMA binaries out (none should be here); include logs and outputs.
            files.append(pth)
    with open(out/cfg["outputs"]["checksums"],"w",encoding="utf-8") as f:
        for pth in files: f.write(f"{sha256_file(pth)}  {pth.relative_to(root)}\n")
    write_text(out/cfg["outputs"]["pass"],f"status=PASS\ntechnical_status=PASS\nstage=MDV6_PGC2019_COMMON_VARIANT_CORRESPONDENCE\nprimary_transfer_n={len(transfer)}\nconditionally_retained_n={len(cond)}\nPGC_role=PRIMARY_TEST_NOT_REPLICATION\nSPARK_inputs_used=NO\nnext=MDV7_SPARK_INDEPENDENT_REPLICATION\n")
    print(summary)
if __name__=="__main__": main()
