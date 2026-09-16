
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

import sys
if len(sys.argv)<2:
    raise SystemExit("Usage: verify_raw26_balance.py /absolute/path/to/frozen_balance_tsv[.gz]")
p=Path(sys.argv[1])
if not p.is_file(): raise SystemExit(f"HOLD: missing {p}")
h,rows=read_tsv(p)
med=pick(h,["median_abs_smd","median_abs_SMD","median_smd","median_absSMD"])
p95=pick(h,["p95_abs_smd","p95_abs_SMD","p95_smd","p95_absSMD"])
if not med or not p95:
    raise SystemExit(f"HOLD: balance file lacks median/P95 SMD columns: {h}")
mvals=[float(z[med]) for z in rows if str(z.get(med,"")).strip()]
pvals=[float(z[p95]) for z in rows if str(z.get(p95,"")).strip()]
if not mvals or not pvals: raise SystemExit("HOLD: no balance values")
status="PASS" if max(mvals)<=0.15 and max(pvals)<=0.25 else "HOLD"
out=[{
 "source_path":str(p),"source_sha256":sha256_file(p),"rows":len(rows),
 "max_median_abs_SMD":max(mvals),"max_p95_abs_SMD":max(pvals),
 "median_gate":"PASS" if max(mvals)<=0.15 else "FAIL",
 "p95_gate":"PASS" if max(pvals)<=0.25 else "FAIL",
 "status":status
}]
write_tsv(OUT/"MAR06_raw26_balance_verified.tsv",list(out[0].keys()),out)
print("MAR06_RAW26_BALANCE="+status)
if status!="PASS": raise SystemExit(40)
