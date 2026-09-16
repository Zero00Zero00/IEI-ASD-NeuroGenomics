
from pathlib import Path
import csv, gzip, json, hashlib, os, math, statistics, subprocess, tarfile, shutil, time
from collections import defaultdict, Counter
from datetime import datetime, timezone

CH3 = Path(os.environ.get("CH3_ROOT","/home/h3021/chapter3"))
MA  = Path(os.environ.get("MA_ROOT","/home/h3021/chapter3/12_molecular_autism_revision"))
PKG = Path(__file__).resolve().parent.parent
OUT = MA / "07_chain_calibration" / "MAR07_chain_calibration_v1.0R4p1"
WP01 = MA / "01_anchor_ingest" / "WP01_raw26_primary_v1"
MAR03 = MA / "03_common_retest" / "MAR03_common_variant_retest_v1.0R2"
MAR04 = MA / "04_effect_boundary" / "MAR04_effect_boundary_v1"
MAR05 = MA / "05_magma_collinearity" / "MAR05_collinearity_v1"
MAR06 = MA / "06_matched_null_audit" / "MAR06_matched_null_audit_v1.0R2"
R3 = MA / "07_chain_calibration" / "MAR07_chain_calibration_v1.0R3"

def utc():
    return datetime.now(timezone.utc).isoformat()

def sha256_file(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()

def open_text(p):
    if str(p).endswith(".gz"):
        return gzip.open(p,"rt",encoding="utf-8-sig",errors="replace")
    return open(p,"rt",encoding="utf-8-sig",errors="replace")

def read_tsv(p):
    with open_text(p) as f:
        r=csv.DictReader(f,delimiter="\t")
        return r.fieldnames or [], list(r)

def write_tsv(p,fields,rows):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True)
    with open(p,"w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

def parse_pass(p):
    d={}
    for line in Path(p).read_text(encoding="utf-8",errors="replace").splitlines():
        if "=" in line:
            k,v=line.split("=",1); d[k]=v
    return d

def wilson_ci(k,n,z=1.959963984540054):
    if n<=0: return (float("nan"),float("nan"))
    p=k/n
    den=1+z*z/n
    center=(p+z*z/(2*n))/den
    half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return max(0.0,center-half),min(1.0,center+half)

def ensure_dirs():
    for d in ["state","logs","observed","certification","calibration","source_data","executor_work"]:
        (OUT/d).mkdir(parents=True,exist_ok=True)

import argparse, pandas as pd, tempfile, yaml, importlib.util
ap=argparse.ArgumentParser()
ap.add_argument("--replicate-id",type=int,required=True); ap.add_argument("--seed",type=int,required=True); ap.add_argument("--out",required=True)
a=ap.parse_args()
cert=OUT/"state"/"MAR07_FULLCHAIN_EXECUTOR_CERTIFIED.txt"
if not cert.is_file(): raise SystemExit("HOLD: executor not certified")
work=OUT/"executor_work"/f"rep_{a.replicate_id:04d}"; work.mkdir(parents=True,exist_ok=True)
start=time.time()
resolved=json.loads((OUT/"MAR07_resolved_inputs.json").read_text())

# Generate pseudo groups.
cp=subprocess.run([os.environ.get("MAR07_PYTHON",os.sys.executable),str(PKG/"executor"/"generate_groups.py"),
                   str(a.replicate_id),str(a.seed),str(work)],
                  stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
(OUT/"logs"/f"rep_{a.replicate_id:04d}_groups.log").write_text(cp.stdout,encoding="utf-8")
if cp.returncode!=0: raise SystemExit("pseudo-group generation failed")

# Create pseudo gene universe for exact MDV2 primary adapter.
g=pd.read_csv(work/"pseudo_groups.tsv",sep="\t")
gu=pd.read_csv(resolved["GENE_UNIVERSE"],sep="\t")
gmap=dict(zip(g["gene"].astype(str),g["R0_group"].astype(str)))
gu["R0_group"]=[gmap[str(x)] for x in gu["HGNC_symbol"].astype(str)]
gu_path=work/"gene_universe.tsv"; gu.to_csv(gu_path,sep="\t",index=False)

mdv2out=work/"mdv2.tsv"
cp=subprocess.run(["Rscript",str(PKG/"certification"/"mdv2_primary_replay.R"),str(gu_path),resolved["PATHWAY_MEMBERSHIP"],resolved["MDV2_CONFIG"],str(mdv2out)],
                  stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
(OUT/"logs"/f"rep_{a.replicate_id:04d}_mdv2.log").write_text(cp.stdout,encoding="utf-8")
if cp.returncode!=0: raise SystemExit("MDV2 failed")

mdv3out=work/"mdv3.tsv"
cp=subprocess.run(["Rscript",str(PKG/"executor"/"mdv3_pseudo.R"),str(CH3),str(work/"pseudo_groups.tsv"),str(work/"mdv3_tmp"),str(mdv3out)],
                  stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
(OUT/"logs"/f"rep_{a.replicate_id:04d}_mdv3.log").write_text(cp.stdout,encoding="utf-8")
if cp.returncode!=0: raise SystemExit("MDV3 failed")

# MDV4 empirical: filter frozen sets so pseudo-overlap genes are absent from own reference null.
groups=pd.read_csv(work/"pseudo_groups.tsv",sep="\t"); pseudo=set(groups.loc[groups.R0_group=="Overlap_raw","gene"].astype(str))
sets=pd.read_csv(resolved["RAW26_SETS"],sep="\t",compression="gzip")
keep=sets.gene_set.map(lambda s: not bool(pseudo.intersection(str(s).split(";"))))
fsets=sets.loc[keep].copy()
if len(fsets) == 0:
    raise SystemExit("HOLD: no pseudo-anchor-independent null sets available")
if len(fsets) < 8000:
    print(
        f"[WARN] reference_B_MDV4={len(fsets)} below provisional "
        "8000 implementation guard; retained under MAR07 Amendment 01"
    )
fset_path=work/"filtered_sets.tsv.gz"; fsets.to_csv(fset_path,sep="\t",index=False,compression="gzip")
tmp=work/"mdv4_tmp"; tmp.mkdir(exist_ok=True)
pd.DataFrame({"gene":sorted(pseudo)}).to_csv(tmp/"WP01_raw26_gene_set.tsv",sep="\t",index=False)
script=Path(resolved["WP01_EMPIRICAL"]); os.sys.path.insert(0,str(script.parent))
spec=importlib.util.spec_from_file_location("wp01_emp_real",script); mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
emp=mod.empirical_for_sets(CH3,tmp,fset_path,"NULL")
md3=pd.read_csv(mdv3out,sep="\t")
x=md3.merge(emp[["source","pathway_id","P_emp","q_emp"]],on=["source","pathway_id"],how="left",validate="one_to_one")
x["Tier1"]=(x["enriched_flag"].astype(str).str.lower().isin(["true","1"]))&(x["q_Firth"]<0.05)&(x["q_emp"]<0.05)
md2=pd.read_csv(mdv2out,sep="\t")
ng2=int(md2["primary_pass"].astype(str).str.lower().isin(["true","1"]).sum())
ng3=int(md3["enriched_flag"].astype(str).str.lower().isin(["true","1"]).sum())
nt=int(x["Tier1"].sum())
bal=pd.read_csv(work/"pseudo_group_balance.tsv",sep="\t")
hashes={}
for grp in ["Overlap_raw","IEI_only","ASD_only"]:
    genes=sorted(groups.loc[groups.R0_group==grp,"gene"].astype(str))
    hashes[grp]=hashlib.sha256(("\n".join(genes)+"\n").encode()).hexdigest()
obj={"replicate_id":a.replicate_id,"seed":a.seed,"n_MDV2":ng2,"n_MDV3":ng3,"n_Tier1":nt,
     "reference_B_MDV4":int(len(fsets)),"group_hashes":hashes,
     "covariate_balance":bal.to_dict("records"),"runtime_sec":time.time()-start}
Path(a.out).write_text(json.dumps(obj,indent=2),encoding="utf-8")
print(json.dumps(obj))
