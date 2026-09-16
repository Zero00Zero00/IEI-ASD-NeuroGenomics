
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

import tarfile
if not (OUT/"MAR05_PASS.txt").is_file(): raise SystemExit("HOLD: MAR05_PASS missing")
handoff=MA/"90_handoff"; handoff.mkdir(parents=True,exist_ok=True)
ts=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
p=handoff/f"MAR05_REVIEW_{ts}.tar.gz"
include=[
 "MAR05_preflight_report.tsv","MAR05_resolved_inputs.json","MAR05_preanalysis_lock.json","MAR05_STOPA_RELEASE.txt",
 "MAR05_membership_reconciliation.tsv","MAR05_unit_overlap_long.tsv","MAR05_design_correlation_long.tsv",
 "MAR05_vif.tsv","MAR05_condition_number.tsv","MAR05_standalone_overlap_long.tsv","MAR05_coefficient_correlation.tsv",
 "MAR05_conditional_model_gate.tsv","MAR05_PGC_jaccard_matrix.tsv","MAR05_SPARK_jaccard_matrix.tsv",
 "MAR05_PGC_phi_matrix.tsv","MAR05_SPARK_phi_matrix.tsv","MAR05_postflight_checks.tsv",
 "MAR05_analysis_lock.json","MAR05_checksums.sha256","MAR05_PASS.txt"
]
with tarfile.open(p,"w:gz") as tf:
    for n in include:
        f=OUT/n
        if f.exists(): tf.add(f,arcname=n)
    for f in (OUT/"source_data").glob("*"): tf.add(f,arcname="source_data/"+f.name)
sha=sha256_file(p)
Path(str(p)+".sha256").write_text(f"{sha}  {p.name}\n",encoding="utf-8")
print(f"REVIEW_PACKAGE={p}")
print(f"SHA256={sha}")
print("PACK_REVIEW=PASS")
