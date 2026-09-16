
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

import sys, importlib.util
ensure_dirs()
resolved=json.loads((OUT/"MAR07_resolved_inputs.json").read_text())
tmp=OUT/"certification"/"mdv4_tmp"
if tmp.exists(): shutil.rmtree(tmp)
tmp.mkdir(parents=True)
shutil.copy2(resolved["RAW26_GENE_SET"],tmp/"WP01_raw26_gene_set.tsv")
script=Path(resolved["WP01_EMPIRICAL"])
sys.path.insert(0,str(script.parent))
spec=importlib.util.spec_from_file_location("wp01_empirical_tiering_real",script)
mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
res=mod.empirical_for_sets(CH3,tmp,Path(resolved["RAW26_SETS"]),"CERT_REPLAY")
_,frozen=read_tsv(Path(resolved["RAW26_TIER"]))
fm={(z["source"],z["pathway_id"]):z for z in frozen}
rm={(str(z["source"]),str(z["pathway_id"])):z for z in res.to_dict("records")}
if set(fm)!=set(rm):
    raise SystemExit(f"HOLD: MDV4 key-universe mismatch frozen={len(fm)} replay={len(rm)}")

def safe_num(v):
    if v is None: return None
    s=str(v).strip()
    if s.upper() in {"","NA","NAN","N/A","NULL","NONE",".","<NA>"}: return None
    try:
        x=float(s)
    except Exception:
        return None
    if math.isnan(x): return None
    return x

rows=[]; fail=0
for k in sorted(fm):
    z=rm[k]; f=fm[k]
    p1=safe_num(z.get("P_emp")); p0=safe_num(f.get("P_emp"))
    q1=safe_num(z.get("q_emp")); q0=safe_num(f.get("q_emp"))
    ok=(p0 is not None and p1 is not None and q0 is not None and q1 is not None)
    dp=abs(p1-p0) if ok else ""
    dq=abs(q1-q0) if ok else ""
    ok=ok and dp<=1e-12 and dq<=1e-12
    status="PASS" if ok else "FAIL"; fail+=int(not ok)
    rows.append({"source":k[0],"pathway_id":k[1],
                 "P_emp_expected":f.get("P_emp",""),"P_emp_replay":z.get("P_emp",""),"delta_P_emp":dp,
                 "q_emp_expected":f.get("q_emp",""),"q_emp_replay":z.get("q_emp",""),"delta_q_emp":dq,
                 "status":status})
write_tsv(OUT/"certification"/"MAR07_MDV4_REPLAY.tsv",list(rows[0].keys()),rows)

tier1=0
for k,z in rm.items():
    f=fm[k]
    orf=safe_num(f.get("OR_Firth"))
    qf=safe_num(f.get("q_Firth"))
    qe=safe_num(z.get("q_emp"))
    is_t1=(orf is not None and qf is not None and qe is not None and orf>1 and qf<0.05 and qe<0.05)
    tier1+=int(is_t1)
if fail or tier1!=39:
    raise SystemExit(f"HOLD: MDV4 replay fail_n={fail}, Tier1={tier1}")
print("MAR07_MDV4_REPLAY=PASS")
