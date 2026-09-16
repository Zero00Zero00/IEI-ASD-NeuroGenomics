
from pathlib import Path
import csv, gzip, json, hashlib, os, math, re, statistics, itertools, shutil, tarfile
from collections import defaultdict, Counter
from datetime import datetime, timezone

MA = Path(os.environ.get("MA_ROOT","/home/h3021/chapter3/12_molecular_autism_revision"))
CH3 = Path(os.environ.get("CH3_ROOT","/home/h3021/chapter3"))
PKG = Path(__file__).resolve().parent.parent
OUT = MA / "06_matched_null_audit" / "MAR06_matched_null_audit_v1.0R2"
WP01 = MA / "01_anchor_ingest" / "WP01_raw26_primary_v1"
MAR03 = MA / "03_common_retest" / "MAR03_common_variant_retest_v1.0R2"
MAR05 = MA / "05_magma_collinearity" / "MAR05_collinearity_v1"
PY = Path(os.environ.get("MAR06_PYTHON", str(MA/"envs"/"mar02_mdv5_legacy"/"bin"/"python")))

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

def iter_tsv(p):
    f=open_text(p)
    r=csv.DictReader(f,delimiter="\t")
    return f,r

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

def pick(header,names):
    m={x.strip().lower():x for x in header}
    for n in names:
        if n.lower() in m: return m[n.lower()]
    return None

def safe_float(x):
    try: return float(x)
    except: return float("nan")

def entropy_effective(counts):
    vals=[float(x) for x in counts if float(x)>0]
    if not vals: return (0.0,0.0,0.0)
    total=sum(vals); p=[x/total for x in vals]
    H=-sum(x*math.log(x) for x in p)
    expH=math.exp(H)
    invsim=1/sum(x*x for x in p)
    return H,expH,invsim

checks=[]
def ck(cid,ok,obs,exp,detail=""):
    checks.append({"check_id":cid,"status":"PASS" if ok else "HOLD","observed":obs,"expected":exp,"detail":detail})

req={
 "WP01_PASS":WP01/"WP01_PASS.txt",
 "MAR03_PASS":MAR03/"MAR03_PASS.txt",
 "MAR05_PASS":MAR05/"MAR05_PASS.txt",
 "MAR03_MATCHED_SETS":MAR03/"matched"/"MAR03_matched_sets.tsv.gz",
 "MAR03_MATCHING_SMD":MAR03/"matched"/"MAR03_matching_smd.tsv.gz",
 "MAR03_CANDIDATE_POOL_QC":MAR03/"matched"/"MAR03_candidate_pool_qc.tsv",
 "MAR03_TILT_DIAG":MAR03/"matched"/"MAR03_matching_tilt_diagnostics.tsv",
 "MAR03_CONTROL_REUSE":MAR03/"matched"/"MAR03_matched_control_reuse.tsv",
 "MAR03_MATCHING_QC":MAR03/"MAR03_matching_qc.tsv",
 "MAR03_MATCHING_SUMMARY":MAR03/"MAR03_matching_summary.tsv",
}
for k,p in req.items(): ck(k,p.is_file(),str(p),"exists")
for k in ["WP01_PASS","MAR03_PASS","MAR05_PASS"]:
    p=req[k]
    if p.is_file(): ck(k+"_status","PASS" in p.read_text(encoding="utf-8",errors="replace"),"PASS marker","PASS")

resf=OUT/"MAR06_RAW26_ARTIFACT_RESOLUTION.json"
ck("RAW26_resolution",resf.is_file(),str(resf),"exists")
rawres={}
if resf.is_file():
    rawres=json.loads(resf.read_text())
    for role in ["primary_K50","sensitivity_K20"]:
        for key,pth in rawres.get(role,{}).items():
            p=Path(pth); ck(f"RAW26_{role}_{key}",p.is_file(),str(p),"exists")

if req["MAR03_MATCHING_SUMMARY"].is_file():
    _,ms=read_tsv(req["MAR03_MATCHING_SUMMARY"])
    ck("common_target_n",len(ms)==20,len(ms),20)
    ck("common_unique_sets_all",all(int(float(z["unique_set_n"]))==10000 for z in ms),
       min(int(float(z["unique_set_n"])) for z in ms),10000)
    ck("common_exact_chr_all",all(z["exact_chromosome_composition"]=="YES" for z in ms),
       sum(z["exact_chromosome_composition"]=="YES" for z in ms),20)
if req["MAR03_MATCHING_QC"].is_file():
    _,mq=read_tsv(req["MAR03_MATCHING_QC"])
    ck("common_QC_all_PASS",all(z["status"]=="PASS" for z in mq),sum(z["status"]=="PASS" for z in mq),len(mq))
ck("validated_python",PY.is_file(),str(PY),"exists")

OUT.mkdir(parents=True,exist_ok=True)
write_tsv(OUT/"MAR06_preflight_report.tsv",["check_id","status","observed","expected","detail"],checks)
resolved={k:str(v) for k,v in req.items()}
resolved.update({"RAW26":rawres,"PYTHON":str(PY)})
(OUT/"MAR06_resolved_inputs.json").write_text(json.dumps(resolved,indent=2),encoding="utf-8")
bad=[z for z in checks if z["status"]!="PASS"]
if bad: raise SystemExit(f"HOLD: MAR06 preflight failed ({len(bad)} checks)")
print("MAR06_PREFLIGHT=PASS")
