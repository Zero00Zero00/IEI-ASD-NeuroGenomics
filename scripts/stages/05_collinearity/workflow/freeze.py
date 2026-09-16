
from pathlib import Path
import csv, json, hashlib, os, math, itertools, shutil, re
from datetime import datetime, timezone

MA = Path(os.environ.get("MA_ROOT","/home/h3021/chapter3/12_molecular_autism_revision"))
CH3 = Path(os.environ.get("CH3_ROOT","/home/h3021/chapter3"))
PKG = Path(__file__).resolve().parent.parent
OUT = MA / "05_magma_collinearity" / "MAR05_collinearity_v1"
MAR03_OUT = MA / "03_common_retest" / "MAR03_common_variant_retest_v1.0R2"
PGC = CH3 / "06_magma_pgc" / "MDV6_v1_1"
SPARK = CH3 / "07_magma_spark" / "MDV7_v1_1"
PY = Path(os.environ.get("MAR05_PYTHON", str(MA/"envs"/"mar02_mdv5_legacy"/"bin"/"python")))

def utc():
    return datetime.now(timezone.utc).isoformat()

def sha256_file(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()

def read_tsv(p):
    with open(p,encoding="utf-8-sig",errors="replace",newline="") as f:
        r=csv.DictReader(f,delimiter="\t")
        return r.fieldnames or [], list(r)

def write_tsv(p,fields,rows):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True)
    with open(p,"w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

def parse_pass(p):
    d={}
    for line in Path(p).read_text(encoding="utf-8").splitlines():
        if "=" in line:
            k,v=line.split("=",1); d[k]=v
    return d

def pick(header,names):
    m={x.strip().lower():x for x in header}
    for n in names:
        if n.lower() in m: return m[n.lower()]
    return None

def resolve_mar04_pass():
    candidates=[
        MA/"MAR04_FREEZE_RELEASE_v1.0R1"/"manifest"/"MAR04_FREEZE_PASS.txt",
        MA/"04_effect_boundary"/"MAR04_effect_boundary_v1"/"MAR04_PASS.txt",
    ]
    good=[p for p in candidates if p.is_file()]
    if not good:
        raise SystemExit("HOLD: MAR04 PASS not found")
    return good[0]

OUT.mkdir(parents=True,exist_ok=True)
pre=OUT/"MAR05_preflight_report.tsv"
if not pre.is_file(): raise SystemExit("HOLD: run preflight first")
_,rows=read_tsv(pre)
if any(z["status"]!="PASS" for z in rows): raise SystemExit("HOLD: preflight not PASS")

lock={
 "stage":"MAR05_MAGMA_COLLINEARITY_AUDIT",
 "version":"1.0R1",
 "timestamp_utc":utc(),
 "primary_family":"A_HIGHER_ORDER_EXPANDED",
 "primary_unit_n":7,
 "secondary_standalone_overlap":"YES",
 "methods":{
   "pairwise_overlap":["intersection_n","union_n","Jaccard","phi_binary_Pearson"],
   "VIF":"For each primary unit indicator, OLS on intercept + other six indicators; VIF=1/(1-R2)",
   "condition_number":"Centered/scaled primary 7-unit binary design matrix; kappa=max(singular)/min(singular)",
   "design_correlation":"Pearson correlation of binary unit indicators (equal to phi)",
   "coefficient_correlation":"Use native MAGMA covariance only if explicitly available; otherwise NOT_AVAILABLE_FROM_NATIVE_OUTPUT",
 },
 "thresholds":{
   "VIF_GREEN_lt":5,
   "VIF_CAUTION_5_to_10":True,
   "VIF_HIGH_gt":10,
   "condition_GREEN_lt":30,
   "condition_CAUTION_30_to_100":True,
   "condition_HIGH_gt":100,
   "pairwise_Jaccard_GREEN_lt":0.5,
   "pairwise_Jaccard_CAUTION_0.5_to_0.7":True,
   "pairwise_Jaccard_HIGH_gt":0.7,
 },
 "conditional_gate":{
   "hard_rule":"If any primary VIF>10 OR condition_number>100 -> SUPPLEMENT_ONLY",
   "otherwise":"MAIN_TEXT_SUPPORTIVE_ELIGIBLE",
   "marginal_model_status":"PRIMARY_ALWAYS"
 },
 "association_coefficients_used":"NO",
 "association_scores_used":"NO",
}
p=OUT/"MAR05_preanalysis_lock.json"
if p.is_file():
    old=json.loads(p.read_text()); old.pop("timestamp_utc",None); new=dict(lock); new.pop("timestamp_utc",None)
    if old!=new: raise SystemExit("HOLD: existing MAR05 lock differs")
    print("MAR05_FREEZE=PASS_REUSED"); raise SystemExit(0)
p.write_text(json.dumps(lock,indent=2),encoding="utf-8")
print("MAR05_FREEZE=PASS")
print("HOLD_FOR_MANUAL_GATE")
