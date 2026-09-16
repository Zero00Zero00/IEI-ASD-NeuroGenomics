
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

import sys, subprocess
checks=[]
def ck(cid,ok,obs,exp,detail=""):
    checks.append({"check_id":cid,"status":"PASS" if ok else "HOLD","observed":obs,"expected":exp,"detail":detail})

req={
 "MAR03_PASS":MAR03_OUT/"MAR03_PASS.txt",
 "MAR03_MEMBERSHIP":MAR03_OUT/"MAR03_frozen_membership.tsv",
 "MAR03_PGC_PRIMARY":MAR03_OUT/"MAR03_pgc_higher_order.tsv",
 "MAR03_SPARK_PRIMARY":MAR03_OUT/"MAR03_spark_higher_order.tsv",
 "PGC_GENE_RESULTS":PGC/"PGC2019_gene_results_w10.tsv",
 "SPARK_GENE_RESULTS":SPARK/"SPARK_EUR_gene_results_w10.tsv",
}
for k,p in req.items(): ck(k,p.is_file(),str(p),"exists")
try:
    mar04p=resolve_mar04_pass()
    ck("MAR04_PASS", "status=PASS" in mar04p.read_text(encoding="utf-8"), str(mar04p),"status=PASS")
except Exception as e:
    mar04p=None; ck("MAR04_PASS",False,repr(e),"status=PASS")

if req["MAR03_PASS"].is_file():
    ck("MAR03_status","status=PASS" in req["MAR03_PASS"].read_text(),"PASS","PASS")

if req["MAR03_MEMBERSHIP"].is_file():
    _,m=read_tsv(req["MAR03_MEMBERSHIP"])
    for ds in ["PGC","SPARK"]:
        for fam,exp in [("A",7),("B",13),("C",7)]:
            n=len({z["target_id"] for z in m if z["dataset"]==ds and z["family"]==fam})
            ck(f"{ds}_{fam}_target_n",n==exp,n,exp)

# Frozen external hashes from the prior real-interface probe.
known={
 str(req["PGC_GENE_RESULTS"]):"f4eb46aca8c3ba205d4d7576e88525e8b2f8c02ec2b9d571dd4282362593106d",
 str(req["SPARK_GENE_RESULTS"]):"99a8bd46f332615262c6e0422be04d3988b154bf6a1af9fa7b5670cf19308702",
}
for p,exp in known.items():
    if Path(p).is_file():
        act=sha256_file(p); ck("hash_"+Path(p).name,act==exp,act,exp)

ck("python_numpy",PY.is_file(),str(PY),"validated Python executable")

OUT.mkdir(parents=True,exist_ok=True)
write_tsv(OUT/"MAR05_preflight_report.tsv",["check_id","status","observed","expected","detail"],checks)
resolved={
 "MAR03_OUT":str(MAR03_OUT),
 "MAR04_PASS":str(mar04p) if mar04p else "",
 "PGC_GENE_RESULTS":str(req["PGC_GENE_RESULTS"]),
 "SPARK_GENE_RESULTS":str(req["SPARK_GENE_RESULTS"]),
 "PYTHON":str(PY),
 "input_hashes":{k:sha256_file(p) for k,p in req.items() if p.is_file()},
}
(OUT/"MAR05_resolved_inputs.json").write_text(json.dumps(resolved,indent=2),encoding="utf-8")
bad=[z for z in checks if z["status"]!="PASS"]
if bad: raise SystemExit(f"HOLD: MAR05 preflight failed ({len(bad)} checks)")
print("MAR05_PREFLIGHT=PASS")
