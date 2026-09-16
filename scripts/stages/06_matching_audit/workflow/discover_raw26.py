
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

authority=json.loads((PKG/"config"/"MAR06_RAW26_AUTHORITY.json").read_text(encoding="utf-8"))
root=Path(authority["wp01_root"])
rows=[]
resolved={}
for role in ["primary_K50","sensitivity_K20"]:
    resolved[role]={}
    for key in ["candidate_pool","assignments","sets","balance_long","balance_summary"]:
        spec=authority[role][key]
        p=root/spec["file"]
        status="PASS" if p.is_file() else "MISSING"
        actual=sha256_file(p) if p.is_file() else ""
        expected=spec.get("sha256") or ""
        hash_status="NOT_FROZEN_IN_PROBE" if not expected else ("PASS" if actual==expected else "FAIL")
        rows.append({"role":role,"artifact":key,"path":str(p),"exists":status,
                     "actual_sha256":actual,"expected_sha256":expected,"hash_status":hash_status})
        if not p.is_file():
            raise SystemExit(f"HOLD: required raw26 artifact missing: {p}")
        if expected and actual!=expected:
            raise SystemExit(f"HOLD: raw26 artifact hash drift: {p}")
        resolved[role][key]=str(p)
resolved["pilot_1k"]={k:str(root/v) for k,v in authority["pilot_1k"].items() if k.endswith("_file")}
write_tsv(OUT/"MAR06_RAW26_AUTHORITY_BINDING.tsv",
          ["role","artifact","path","exists","actual_sha256","expected_sha256","hash_status"],rows)
(OUT/"MAR06_RAW26_ARTIFACT_RESOLUTION.json").write_text(json.dumps(resolved,indent=2),encoding="utf-8")
print("MAR06_RAW26_AUTHORITY=PASS")
print("PRIMARY_K50_POOL=",resolved["primary_K50"]["candidate_pool"])
print("PRIMARY_K50_ASSIGNMENTS=",resolved["primary_K50"]["assignments"])
print("PRIMARY_K50_SETS=",resolved["primary_K50"]["sets"])
print("PRIMARY_K50_BALANCE=",resolved["primary_K50"]["balance_long"])
print("K20_SENSITIVITY=BOUND")
