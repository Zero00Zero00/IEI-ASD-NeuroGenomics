
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

checks=[]
def ck(cid,ok,obs,exp,detail=""):
    checks.append({"check_id":cid,"status":"PASS" if ok else "HOLD","observed":obs,"expected":exp,"detail":detail})
req=["MAR05_membership_reconciliation.tsv","MAR05_unit_overlap_long.tsv","MAR05_design_correlation_long.tsv",
     "MAR05_vif.tsv","MAR05_condition_number.tsv","MAR05_standalone_overlap_long.tsv",
     "MAR05_coefficient_correlation.tsv","MAR05_conditional_model_gate.tsv",
     "MAR05_PGC_jaccard_matrix.tsv","MAR05_SPARK_jaccard_matrix.tsv",
     "MAR05_PGC_phi_matrix.tsv","MAR05_SPARK_phi_matrix.tsv"]
for f in req: ck("exists_"+f,(OUT/f).is_file(),f,"exists")
if all((OUT/f).is_file() for f in req):
    _,rec=read_tsv(OUT/"MAR05_membership_reconciliation.tsv")
    _,pair=read_tsv(OUT/"MAR05_unit_overlap_long.tsv")
    _,vif=read_tsv(OUT/"MAR05_vif.tsv")
    _,cond=read_tsv(OUT/"MAR05_condition_number.tsv")
    _,spair=read_tsv(OUT/"MAR05_standalone_overlap_long.tsv")
    ck("membership_reconcile_n",len(rec)==14,len(rec),14)
    ck("membership_reconcile_all",all(z["match"]=="YES" for z in rec),sum(z["match"]=="YES" for z in rec),14)
    ck("primary_pair_n",len(pair)==42,len(pair),42)
    ck("vif_n",len(vif)==14,len(vif),14)
    ck("condition_n",len(cond)==2,len(cond),2)
    ck("standalone_pair_n",len(spair)==156,len(spair),156)
    ck("vif_finite",all(math.isfinite(float(z["VIF"])) for z in vif),"","finite")
    ck("condition_finite",all(math.isfinite(float(z["condition_number"])) for z in cond),"","finite")

    # Frozen primary expectation oracle; compare only if input authorities match.
    _,ref=read_tsv(PKG/"config"/"MAR05_PRIMARY_REFERENCE_EXPECTATION.tsv")
    # pair/VIF/condition rounded tolerance
    tol=1e-8; ok=True; details=[]
    for ds in ["PGC","SPARK"]:
        actual_pairs={(z["unit_i"],z["unit_j"]):z for z in pair if z["dataset"]==ds}
        for r in [x for x in ref if x["dataset"]==ds and x["metric_type"]=="PAIR"]:
            a=actual_pairs.get((r["unit_i"],r["unit_j"]))
            if not a or abs(float(a["Jaccard"])-float(r["Jaccard"]))>tol or abs(float(a["phi"])-float(r["phi"]))>tol:
                ok=False; details.append(f"{ds}:{r['unit_i']}-{r['unit_j']}")
        actual_v={z["unit_id"]:z for z in vif if z["dataset"]==ds}
        for r in [x for x in ref if x["dataset"]==ds and x["metric_type"]=="VIF"]:
            a=actual_v.get(r["unit_id"])
            if not a or abs(float(a["VIF"])-float(r["VIF"]))>tol:
                ok=False; details.append(f"{ds}:VIF:{r['unit_id']}")
        ar=next(z for z in cond if z["dataset"]==ds)
        rr=next(x for x in ref if x["dataset"]==ds and x["metric_type"]=="CONDITION")
        if abs(float(ar["condition_number"])-float(rr["condition_number"]))>tol:
            ok=False; details.append(f"{ds}:condition")
    ck("primary_reference_replay",ok,";".join(details) if details else "exact within 1e-8","PASS")

write_tsv(OUT/"MAR05_postflight_checks.tsv",["check_id","status","observed","expected","detail"],checks)
bad=[z for z in checks if z["status"]!="PASS"]
if bad: raise SystemExit(f"HOLD: MAR05 postflight {len(bad)} failures")

# Canonical source-data.
sd=OUT/"source_data"; sd.mkdir(exist_ok=True)
for s,d in [
 ("MAR05_unit_overlap_long.tsv","MAR05_SOURCE_UNIT_OVERLAP_LONG.tsv"),
 ("MAR05_vif.tsv","MAR05_SOURCE_VIF.tsv"),
 ("MAR05_condition_number.tsv","MAR05_SOURCE_CONDITION_NUMBER.tsv"),
 ("MAR05_design_correlation_long.tsv","MAR05_SOURCE_DESIGN_CORRELATION.tsv"),
 ("MAR05_standalone_overlap_long.tsv","MAR05_SOURCE_STANDALONE_OVERLAP.tsv"),
 ("MAR05_coefficient_correlation.tsv","MAR05_SOURCE_COEFFICIENT_CORRELATION.tsv"),
 ("MAR05_conditional_model_gate.tsv","MAR05_SOURCE_CONDITIONAL_MODEL_GATE.tsv"),
 ("MAR05_PGC_jaccard_matrix.tsv","MAR05_PLOT_PGC_JACCARD_MATRIX.tsv"),
 ("MAR05_SPARK_jaccard_matrix.tsv","MAR05_PLOT_SPARK_JACCARD_MATRIX.tsv"),
 ("MAR05_PGC_phi_matrix.tsv","MAR05_PLOT_PGC_PHI_MATRIX.tsv"),
 ("MAR05_SPARK_phi_matrix.tsv","MAR05_PLOT_SPARK_PHI_MATRIX.tsv"),
]:
    shutil.copy2(OUT/s,sd/d)

_,gate=read_tsv(OUT/"MAR05_conditional_model_gate.tsv")
combined=next(z for z in gate if z["dataset"]=="COMBINED")
summary=[
 {"metric":"MAR05_status","value":"PASS","note":"Technical audit complete"},
 {"metric":"max_VIF","value":combined["max_VIF"],"note":"Hard threshold >10"},
 {"metric":"max_condition_number","value":combined["condition_number"],"note":"Hard threshold >100"},
 {"metric":"max_pairwise_Jaccard","value":combined["max_pairwise_Jaccard"],"note":"High-overlap diagnostic >0.7"},
 {"metric":"hard_high_collinearity","value":combined["hard_high_collinearity"],"note":""},
 {"metric":"marginal_model_status","value":"PRIMARY","note":"Always primary"},
 {"metric":"conditional_reporting_gate","value":combined["conditional_reporting_gate"],"note":"Main-text supportive vs supplement-only"},
 {"metric":"next_stage","value":"MAR06","note":""},
]
write_tsv(sd/"MAR05_SOURCE_WRITEUP_SUMMARY.tsv",list(summary[0].keys()),summary)

# Locks/checksums/pass.
files=[p for p in OUT.rglob("*") if p.is_file() and p.name not in {"MAR05_checksums.sha256","MAR05_PASS.txt"}]
lock={"stage":"MAR05","status":"PASS","timestamp_utc":utc(),
      "preanalysis_lock_sha256":sha256_file(OUT/"MAR05_preanalysis_lock.json"),
      "output_hashes":{str(p.relative_to(OUT)):sha256_file(p) for p in files}}
(OUT/"MAR05_analysis_lock.json").write_text(json.dumps(lock,indent=2),encoding="utf-8")
lines=[f"{sha256_file(p)}  {p.relative_to(MA)}" for p in sorted([p for p in OUT.rglob("*") if p.is_file() and p.name!="MAR05_checksums.sha256"])]
(OUT/"MAR05_checksums.sha256").write_text("\n".join(lines)+"\n",encoding="utf-8")
(OUT/"MAR05_PASS.txt").write_text(
    "status=PASS\nstage=MAR05_MAGMA_COLLINEARITY_AUDIT\n"
    f"hard_high_collinearity={combined['hard_high_collinearity']}\n"
    "marginal_model_status=PRIMARY\n"
    f"conditional_reporting_gate={combined['conditional_reporting_gate']}\n"
    "next_stage=MAR06\n",encoding="utf-8"
)
print("MAR05_POSTFLIGHT=PASS")
