#!/usr/bin/env python3
import argparse, json, platform, sys
from pathlib import Path
import numpy as np, pandas as pd, scipy, yaml
from mdv7_common import *

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); a=ap.parse_args(); cfg=load_cfg(a.config); root=Path(cfg["project_root"]); out=out_dir(cfg); prov=prov_dir(cfg); checks=[]
    def ck(name,ok,obs="",exp=""): checks.append({"check":name,"status":"PASS" if ok else "FAIL","observed":obs,"expected":exp})
    for k in ["preflight_pass","preanalysis_pass","harmonization_pass","primary_gene_pass","expanded_pass","coreseed_empirical_pass","compact_pass","sensitivity_pass","concordance_pass","precision_pass"]:
        p=out/cfg["outputs"][k]; ck(k,p.exists() and p.stat().st_size>0,str(p),"existing non-empty")
    for gate in ["MDV7_STOPA_RELEASE.txt","MDV7_STOPB_RELEASE.txt"]: ck(gate,(prov/gate).exists() and (prov/gate).stat().st_size>0,str(prov/gate),"manual release exists")
    for st in ["mdv0","mdv1","mdv2","mdv3","mdv4","mdv5","mdv6"]:
        vr=verify_sha256_manifest(root_path(cfg,cfg["upstream"][f"{st}_checksums"]),root); bad=sum(x["status"]!="PASS" for x in vr); ck(f"{st}_checksums_terminal",bad==0,f"failed={bad}",0)
    lock=json.load(open(out/cfg["outputs"]["preanalysis_lock"]))
    for rel,h in lock["scientific_files"].items():
        fp=root/rel; got=sha256_file(fp) if fp.exists() else "MISSING"; ck("freeze_"+Path(rel).name,got==h,got,h)
    pkg=Path(a.config).parent.parent
    for rel,h in lock["code_and_config_sha256"].items():
        p=pkg/rel; got=sha256_file(p) if p.exists() else "MISSING"; ck("code_"+rel,got==h,got,h)
    p=pd.read_csv(out/cfg["outputs"]["expanded_summary"],sep="\t"); ck("primary_unit_n",len(p)==7,len(p),7); ck("primary_ids",sorted(p.unit_id.tolist())==sorted(cfg["expected"]["unit_ids"]),";".join(sorted(p.unit_id)),";".join(sorted(cfg["expected"]["unit_ids"]))); ck("primary_BH",np.allclose(bh(p.P_marginal),p.q_marginal,atol=1e-12,rtol=0)); ck("conditional_BH",np.allclose(bh(p.P_conditional),p.q_conditional,atol=1e-12,rtol=0))
    c=pd.read_csv(out/cfg["outputs"]["compact_summary"],sep="\t"); ck("compact_unit_n",len(c)==7,len(c),7); ck("compact_ids",sorted(c.unit_id.tolist())==sorted(cfg["expected"]["unit_ids"])); ck("compact_BH",np.allclose(bh(c.P_marginal),c.q_marginal,atol=1e-12,rtol=0))
    res=json.load(open(out/cfg["outputs"]["resolved_inputs"])); s=json.dumps(res).lower(); ck("exact_SPARK_only_EUR",Path(res["spark_primary"]).name=="SPARK_update_model1QC_EUR_only.tsv.gz",Path(res["spark_primary"]).name,"SPARK_update_model1QC_EUR_only.tsv.gz"); ck("forbidden_meta_absent",all(Path(cfg["spark"][k]).name.lower() not in s for k in ["excluded_meta_eur","excluded_meta_all"])); ck("multiancestry_not_inferential",Path(cfg["spark"]["multiancestry_sensitivity_candidate"]).name.lower() not in s)
    h=pd.read_csv(out/cfg["outputs"]["liftover_audit"],sep="\t").iloc[0]; ck("liftover_success_gate",float(h.liftover_success_rate)>=cfg["harmonization"]["liftover_success_min"],h.liftover_success_rate,cfg["harmonization"]["liftover_success_min"]); ck("retention_catastrophic_gate",float(h.reference_compatible_retention)>=cfg["harmonization"]["reference_compatible_retention_min"],h.reference_compatible_retention,cfg["harmonization"]["reference_compatible_retention_min"])
    ss=pd.read_csv(out/cfg["outputs"]["sensitivity_summary"],sep="\t"); ck("sensitivity_scenarios",set(ss.scenario)=={"compact","w0","w50","mhc","remove_coreseed","alltier1","nested"},";".join(sorted(set(ss.scenario))))
    cc=pd.read_csv(out/cfg["outputs"]["cross_cohort"],sep="\t"); ck("concordance_unit_n",len(cc)==7,len(cc),7); pr=pd.read_csv(out/cfg["outputs"]["precision_audit"],sep="\t"); ck("precision_unit_n",len(pr)==7,len(pr),7)
    cdf=pd.DataFrame(checks); cdf.to_csv(out/cfg["outputs"]["postflight_checks"],sep="\t",index=False); failed=int((cdf.status=="FAIL").sum())
    write_text(out/cfg["outputs"]["session_versions"],f"python={sys.version}\npandas={pd.__version__}\nnumpy={np.__version__}\nscipy={scipy.__version__}\npyyaml={yaml.__version__}\nplatform={platform.platform()}\nmagma={res.get('magma_version','')}\n")
    if failed: raise RuntimeError(f"MDV7 postflight failed {failed} checks")
    nprim=int(pd.to_numeric(p.primary_support_flag,errors="coerce").fillna(0).astype(bool).sum()); ncond=int(pd.to_numeric(p.conditional_support_flag,errors="coerce").fillna(0).astype(bool).sum()); m02=c[c.unit_id=="M02"].iloc[0]; m02f=bool(m02.secondary_FDR_flag)
    summary=("MDV-7 SPARK INDEPENDENT COMMON-VARIANT FOLLOW-UP — VALIDATION SUMMARY\n"+"="*76+"\n"+f"version={cfg['version']}\nprimary_dataset=SPARK_ONLY_EUR\nprimary_window_kb=10\nprimary_gene_sets=Domain-Expanded\nprimary_unit_n=7\nSPARK_primary_support_n={nprim}\nconditionally_supported_n={ncond}\nM02_Compact_secondary_FDR_positive={m02f}\nreplication_language_allowed=NO\nmultiancestry_inferential_use=NO\nforbidden_SPARK_PGC_meta_files_used=NO\ntechnical_checks_passed_n={len(cdf)}\ntechnical_checks_failed_n=0\n")
    write_text(out/cfg["outputs"]["validation_summary"],summary)
    alock={"version":cfg["version"],"stage":cfg["stage"],"association_rows_read":"YES","SPARK_role":"INDEPENDENT_COMMON_VARIANT_FOLLOWUP_NOT_REPLICATION","primary_dataset":"SPARK_ONLY_EUR","primary_window_kb":10,"primary_set":"Domain-Expanded","SPARK_primary_support_units":p.loc[p.primary_support_flag.astype(bool),"unit_id"].tolist(),"SPARK_conditionally_supported_units":p.loc[p.conditional_support_flag.astype(bool),"unit_id"].tolist(),"M02_Compact_secondary_FDR_positive":m02f,"input_lock_sha256":sha256_file(out/cfg["outputs"]["input_lock"]),"preanalysis_lock_sha256":sha256_file(out/cfg["outputs"]["preanalysis_lock"]),"STOPA_release_sha256":sha256_file(prov/"MDV7_STOPA_RELEASE.txt"),"STOPB_release_sha256":sha256_file(prov/"MDV7_STOPB_RELEASE.txt"),"interpretation_lock":cfg["interpretation_lock"]}; write_json(out/cfg["outputs"]["analysis_lock"],alock)
    exclude={cfg["outputs"]["checksums"],cfg["outputs"]["pass"]}; files=[x for x in sorted(out.rglob("*")) if x.is_file() and x.name not in exclude]
    with open(out/cfg["outputs"]["checksums"],"w") as f:
        for x in files: f.write(f"{sha256_file(x)}  {x.relative_to(root)}\n")
    write_text(out/cfg["outputs"]["pass"],f"status=PASS\ntechnical_status=PASS\nstage=MDV7_SPARK_INDEPENDENT_COMMON_VARIANT_FOLLOWUP\nSPARK_primary_support_n={nprim}\nconditionally_supported_n={ncond}\nM02_Compact_secondary_FDR_positive={m02f}\nSPARK_role=INDEPENDENT_FOLLOWUP_NOT_REPLICATION\nnext=MANUSCRIPT_INTEGRATION_OR_PRESPECIFIED_MDV8_ONLY\n")
    print(summary)
if __name__=="__main__": main()
