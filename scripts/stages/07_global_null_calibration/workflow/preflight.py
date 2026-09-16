
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

ensure_dirs()
checks=[]
def ck(cid,ok,obs,exp,detail=""):
    checks.append({"check_id":cid,"status":"PASS" if ok else "HOLD","observed":obs,"expected":exp,"detail":detail})

req={
 "WP01_PASS":WP01/"WP01_PASS.txt",
 "MODEL_FRAME":WP01/"WP01_model_frame.tsv.gz",
 "RAW26_GENE_SET":WP01/"WP01_raw26_gene_set.tsv",
 "RAW26_MDV3":WP01/"WP01_raw26_mdv3_all.tsv.gz",
 "RAW26_TIER":WP01/"WP01_raw26_tier.tsv.gz",
 "RAW26_TIER1":WP01/"WP01_raw26_tier1.tsv",
 "RAW26_SETS":WP01/"WP01_raw26_full_sets_10k.tsv.gz",
 "PATHWAY_MEMBERSHIP":CH3/"02_pathways"/"03_pathway_membership.tsv.gz",
 "GENE_UNIVERSE":CH3/"01_gene_universe"/"01_gene_universe.tsv",
 "GENE_GROUPS":CH3/"01_gene_universe"/"02_gene_groups.tsv",
 "MDV2_FROZEN":CH3/"03_intersection"/"03_intersection_results.tsv.gz",
 "MDV2_MAIN":CH3/"workflow"/"scripts"/"mdv2_01_dual_firth_main.R",
 "MDV2_CONFIG":CH3/"config"/"mdv2_dual_firth_gate_v1.yaml",
 "WP01_MDV3_BRIDGE":MA/"workflow"/"scripts"/"wp01_mdv3_legacy_bridge.R",
 "LEGACY_SOLVER":CH3/"workflow"/"scripts"/"mdv3_firth_solver_v1p3.R",
 "MDV3_CONFIG":CH3/"config"/"mdv3_coreseed_bias_firth_gate_v1.yaml",
 "WP01_EMPIRICAL":MA/"workflow"/"scripts"/"wp01_empirical_tiering.py",
 "WP01_COMMON":MA/"workflow"/"scripts"/"wp01_common.py",
 "WP01_GLOBAL":MA/"workflow"/"scripts"/"wp01_global_overlap.py",
 "MAR03_PASS":MAR03/"MAR03_PASS.txt",
 "MAR04_PASS":MAR04/"MAR04_PASS.txt",
 "MAR05_PASS":MAR05/"MAR05_PASS.txt",
 "MAR06_PASS":MAR06/"MAR06_PASS.txt",
}
for k,p in req.items(): ck(k,p.is_file(),str(p),"exists")
for k in ["WP01_PASS","MAR03_PASS","MAR04_PASS","MAR05_PASS","MAR06_PASS"]:
    p=req[k]
    if p.is_file(): ck(k+"_status","PASS" in p.read_text(encoding="utf-8",errors="replace"),"PASS marker","PASS")

known_hashes={
 str(req["MODEL_FRAME"]):"544b05e8b14e480e1015066480c8064ac80a23c0fa248a3ead034344643b49dc",
 str(req["RAW26_SETS"]):"aeddbb1a841d60475301993d770a973770403b082027d3fc2c558d76c188301a",
 str(req["PATHWAY_MEMBERSHIP"]):"fbf5034608b0397647652fea21816abbe4eaa214dd800a7fdaf9ff4660caa3da",
 str(req["GENE_UNIVERSE"]):"451068856fd057895729e32d64fdc29c098f095b55d73e0d2e1454323f11fce4",
 str(req["MDV2_MAIN"]):"2f9cef3bcfbd8689cc9a67333fdb5056a1435fe69e1acd580886c0538937bfd9",
 str(req["MDV2_CONFIG"]):"f8ea243bf593bde4502ba7e65add2e5e47962c4356c419edcb7b269ca29f9a3a",
 str(req["WP01_MDV3_BRIDGE"]):"7cd9ff6796ec9e513bdadc1ca0aee9184f72335a217b04ada58976eda31e1abb",
 str(req["LEGACY_SOLVER"]):"6f41706188865f6c88e5774855604a62d54fa948c39a8f7fd3032d8efa55d044",
 str(req["MDV3_CONFIG"]):"68f575d01f0149c613bf2c66b691e6823a27275f8ba63471d71750a7965beaa6",
 str(req["WP01_EMPIRICAL"]):"170a984aa1f339f334f5ef7972a552a0b6e4a93057c3aa26872f2c71b28bed54",
 str(req["WP01_COMMON"]):"4e14d150312a6aba382e90ea13a077ce242ff924eeb15693392edf03c28f25f2",
 str(req["WP01_GLOBAL"]):"42ba8c2dc34341925ec3c1a09b57dec1f72cf0ce5bfcf285c28aebf46a804b0c",
}
for p,exp in known_hashes.items():
    pp=Path(p)
    if pp.is_file():
        act=sha256_file(pp); ck("hash_"+pp.name,act==exp,act,exp)

# Runtime contract for certified MDV3.
try:
    cp=subprocess.run(["Rscript","-e",
        'cat(R.version.string,"\\n"); cat(as.character(packageVersion("brglm2")),"\\n")'],
        stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=30)
    lines=[x.strip() for x in cp.stdout.splitlines() if x.strip()]
    ck("R_4_5_1",any("R version 4.5.1" in x for x in lines),lines,"R 4.5.1")
    ck("brglm2_1_1_0",any(x=="1.1.0" for x in lines),lines,"1.1.0")
except Exception as e:
    ck("R_runtime",False,repr(e),"R 4.5.1 + brglm2 1.1.0")

# Observed R3 availability.
for name in [
 "observed/MAR07_bh_by_sensitivity.tsv",
 "observed/MAR07_observed_tier1_stability.tsv",
 "observed/MAR07_observed_dependency_summary.tsv"]:
    p=R3/name; ck("R3_"+Path(name).name,p.is_file(),str(p),"exists")

write_tsv(OUT/"MAR07_preflight_report.tsv",["check_id","status","observed","expected","detail"],checks)
resolved={k:str(v) for k,v in req.items()}
resolved["hashes"]={k:sha256_file(v) for k,v in req.items() if v.is_file()}
(OUT/"MAR07_resolved_inputs.json").write_text(json.dumps(resolved,indent=2),encoding="utf-8")
bad=[z for z in checks if z["status"]!="PASS"]
if bad: raise SystemExit(f"HOLD: MAR07 R4 preflight failed ({len(bad)} checks)")
print("MAR07_R4_PREFLIGHT=PASS")
