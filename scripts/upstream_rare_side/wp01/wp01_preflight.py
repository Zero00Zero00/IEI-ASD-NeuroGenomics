#!/usr/bin/env python3
from pathlib import Path
import json, os, subprocess, sys, pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from wp01_common import *

def main():
    if len(sys.argv)<2:
        print("usage: wp01_preflight.py CH3_ROOT",file=sys.stderr); return 2
    root=Path(sys.argv[1]).resolve()
    ma=root/"12_molecular_autism_revision"
    out=ma/"01_anchor_ingest/WP01_raw26_primary_v1"
    out.mkdir(parents=True,exist_ok=True)
    rows=[]
    def add(check,status,critical,details=""):
        rows.append({"check_id":check,"status":status,"critical":"YES" if critical else "NO","details":details})
    # R00 closure
    r00=ma/"00_contract/MAR00_PASS.txt"
    add("R00_PASS", "PASS" if r00.exists() and "PASS" in r00.read_text() else "FAIL", True, str(r00))
    val=ma/"00_contract/MAR00_validation_summary.txt"
    vtxt=val.read_text(errors="replace") if val.exists() else ""
    add("R00_validation_status", "PASS" if "status=PASS" in vtxt and "legacy_changed=0" in vtxt else "FAIL", True, vtxt[:300].replace("\n"," | "))
    sem=ma/"00_contract/MAR00_semantic_metrics.tsv"
    if sem.exists():
        sdf=pd.read_csv(sem,sep="\t")
        bad=sdf[(sdf["blocker"]=="YES") & (sdf["status"]!="PASS")]
        add("R00_semantic_blockers","PASS" if len(bad)==0 else "FAIL",True,f"nonpass_blockers={len(bad)}")
    else: add("R00_semantic_blockers","FAIL",True,"missing semantic metrics")

    files=[
      ("Gene_universe","01_gene_universe/01_gene_universe.tsv"),
      ("Gene_groups","01_gene_universe/02_gene_groups.tsv"),
      ("Pathway_manifest","02_pathways/pathway_manifest.tsv"),
      ("Pathway_membership","02_pathways/03_pathway_membership.tsv.gz"),
      ("Annotation_degree","02_pathways/annotation_degree.tsv"),
      ("Legacy_MDV3_all","04_bias_empirical/04_pathway_all.tsv.gz"),
      ("Legacy_MDV3_spline","04_bias_empirical/MDV3_spline_basis.tsv.gz"),
      ("Legacy_MDV4_tier_all","04_bias_empirical/MDV4_empirical/MDV4_tier_all.tsv.gz"),
      ("Legacy_MDV4_tier1","04_bias_empirical/MDV4_empirical/05_tier1_pathways.tsv"),
      ("Legacy_MDV3_PASS","04_bias_empirical/MDV3_PASS.txt"),
      ("Legacy_MDV4_PASS","04_bias_empirical/MDV4_empirical/MDV4_PASS.txt"),
    ]
    for role,rel in files:
        p=root/rel
        add("INPUT::"+role,"PASS" if p.exists() and p.stat().st_size>0 else "FAIL",True,f"{p}; size={p.stat().st_size if p.exists() else 'NA'}")
    # schema
    try:
        frame,meta=load_structural_frame(root)
        add("SCHEMA::structural_frame","PASS",True,f"rows={len(frame)}; {meta}")
        add("COUNT::universe","PASS" if len(frame)==19267 else "FAIL",True,f"rows={len(frame)}")
        counts={
          "iuis":int(frame.IUIS_flag.sum()),"sfari":int(frame.SFARI_R0_highconf_flag.sum()),
          "raw26":int((frame.IUIS_flag & frame.SFARI_R0_highconf_flag).sum()),
          "core25":int(frame.CoreSeed_flag.sum())
        }
        add("COUNT::gene_groups","PASS" if counts=={"iuis":500,"sfari":935,"raw26":26,"core25":25} else "FAIL",True,str(counts))
    except Exception as e:
        add("SCHEMA::structural_frame","FAIL",True,repr(e))
    try:
        pm=load_pathway_manifest(root); mem=load_pathway_membership(root)
        nsrc=pm.groupby("source").size().to_dict()
        add("SCHEMA::pathways","PASS" if len(pm)==6671 else "FAIL",True,f"manifest_primary={len(pm)}; sources={nsrc}; membership_paths={mem[['source','pathway_id']].drop_duplicates().shape[0]}")
    except Exception as e:
        add("SCHEMA::pathways","FAIL",True,repr(e))
    # Runtime
    for cmd in ["python3","Rscript","bash","tar","sha256sum"]:
        p=subprocess.run(["bash","-lc",f"command -v {cmd}"],capture_output=True,text=True)
        add("RUNTIME::"+cmd,"PASS" if p.returncode==0 else "FAIL",True,p.stdout.strip() or p.stderr.strip())
    try:
        import numpy, pandas
        add("PYTHON::numpy_pandas","PASS",True,f"numpy={numpy.__version__}; pandas={pandas.__version__}")
    except Exception as e: add("PYTHON::numpy_pandas","FAIL",True,repr(e))
    rcheck=subprocess.run(["bash","-lc","Rscript -e 'pkgs<-c(\"data.table\",\"brglm2\",\"mgcv\"); x<-sapply(pkgs,requireNamespace,quietly=TRUE); cat(paste(names(x),x,sep=\"=\",collapse=\";\")); quit(status=if(all(x))0 else 2)'"],capture_output=True,text=True)
    add("R_PACKAGES::data.table_brglm2_mgcv","PASS" if rcheck.returncode==0 else "FAIL",True,(rcheck.stdout+rcheck.stderr).strip())
    df=pd.DataFrame(rows)
    df.to_csv(out/"WP01_preflight_report.tsv",sep="\t",index=False)
    crit=df[(df.critical=="YES") & (df.status!="PASS")]
    summary={"status":"PASS" if len(crit)==0 else "FAIL","critical_failures":len(crit),"checks_total":len(df)}
    write_json(out/"WP01_preflight_summary.json",summary)
    marker=out/("WP01_PREFLIGHT_PASS.txt" if len(crit)==0 else "WP01_PREFLIGHT_HOLD.txt")
    marker.write_text(("PASS" if len(crit)==0 else "HOLD")+"\n",encoding="utf-8")
    print(json.dumps(summary,indent=2))
    return 0 if len(crit)==0 else 3
if __name__=="__main__": raise SystemExit(main())
