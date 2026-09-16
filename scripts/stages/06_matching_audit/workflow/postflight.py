
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

req=["MAR06_raw26_candidate_support.tsv","MAR06_raw26_reuse_summary.tsv","MAR06_raw26_set_integrity.tsv",
     "MAR06_raw26_balance.tsv","MAR06_common_matching_qc.tsv","MAR06_common_control_reuse.tsv","MAR06_common_candidate_pool_qc.tsv"]
for f in req: ck("exists_"+f,(OUT/f).is_file(),f,"exists")
if all((OUT/f).is_file() for f in req):
    _,reuse=read_tsv(OUT/"MAR06_raw26_reuse_summary.tsv")
    _,integ=read_tsv(OUT/"MAR06_raw26_set_integrity.tsv")
    _,bal=read_tsv(OUT/"MAR06_raw26_balance.tsv")
    _,cq=read_tsv(OUT/"MAR06_common_matching_qc.tsv")
    for label,K in [("PRIMARY_K50",50),("SENSITIVITY_K20",20)]:
        rr=[z for z in reuse if z["analysis"]==label]
        ii=next(z for z in integ if z["analysis"]==label)
        bb=[z for z in bal if z["analysis"]==label]
        ck(label+"_anchor_n",len(rr)==26,len(rr),26)
        ck(label+"_candidate_pool",int(ii["candidate_pool_per_anchor_min"])==K and int(ii["candidate_pool_per_anchor_max"])==K,
           f"{ii['candidate_pool_per_anchor_min']}/{ii['candidate_pool_per_anchor_max']}",K)
        ck(label+"_B",int(ii["accepted_replicates"])==10000,ii["accepted_replicates"],10000)
        ck(label+"_unique_sets",int(ii["unique_set_hashes"])==10000,ii["unique_set_hashes"],10000)
        ck(label+"_within_set_duplicates",int(ii["within_set_duplicate_sets"])==0,ii["within_set_duplicate_sets"],0)
        ck(label+"_assignment_membership",int(ii["assignment_set_membership_mismatch"])==0,ii["assignment_set_membership_mismatch"],0)
        ck(label+"_all_candidates_used",all(int(z["selected_unique_candidates"])==K for z in rr),
           min(int(z["selected_unique_candidates"]) for z in rr),K)
        ck(label+"_balance",all(z["status"]=="PASS" for z in bb),
           f"{sum(z['status']=='PASS' for z in bb)}/{len(bb)}","all PASS")
    ck("common_target_n",len(cq)==20,len(cq),20)
    ck("common_B",all(int(float(z["accepted_sets"]))==10000 for z in cq),min(int(float(z["accepted_sets"])) for z in cq),10000)
    ck("common_unique_sets",all(int(float(z["unique_sets"]))==10000 for z in cq),min(int(float(z["unique_sets"])) for z in cq),10000)
    ck("common_exact_chr",all(z["exact_chromosome_composition"]=="YES" for z in cq),sum(z["exact_chromosome_composition"]=="YES" for z in cq),20)
    ck("common_balance_median",all(float(z["max_median_abs_SMD"])<=0.15 for z in cq if z["max_median_abs_SMD"]!=""),
       max(float(z["max_median_abs_SMD"]) for z in cq if z["max_median_abs_SMD"]!=""),"<=0.15")
    ck("common_balance_p95",all(float(z["max_p95_abs_SMD"])<=0.25 for z in cq if z["max_p95_abs_SMD"]!=""),
       max(float(z["max_p95_abs_SMD"]) for z in cq if z["max_p95_abs_SMD"]!=""),"<=0.25")

write_tsv(OUT/"MAR06_postflight_checks.tsv",["check_id","status","observed","expected","detail"],checks)
bad=[z for z in checks if z["status"]!="PASS"]
if bad:
    print(f"HOLD: MAR06 postflight {len(bad)} failure(s); see MAR06_postflight_checks.tsv")
    raise SystemExit(40)

_,reuse=read_tsv(OUT/"MAR06_raw26_reuse_summary.tsv")
_,cq=read_tsv(OUT/"MAR06_common_matching_qc.tsv")
primary=[z for z in reuse if z["analysis"]=="PRIMARY_K50"]
k20=[z for z in reuse if z["analysis"]=="SENSITIVITY_K20"]
gate=[
 {"component":"RAW26_PRIMARY_K50","hard_gate":"PASS",
  "max_reuse_prop":max(float(z["max_reuse_prop"]) for z in primary),
  "minimum_inverse_simpson_effective_support":min(float(z["effective_support_inverse_simpson"]) for z in primary),
  "reuse_status":"DIAGNOSTIC"},
 {"component":"RAW26_SENSITIVITY_K20","hard_gate":"PASS",
  "max_reuse_prop":max(float(z["max_reuse_prop"]) for z in k20),
  "minimum_inverse_simpson_effective_support":min(float(z["effective_support_inverse_simpson"]) for z in k20),
  "reuse_status":"SENSITIVITY_DIAGNOSTIC"},
 {"component":"COMMON_SIDE","hard_gate":"PASS",
  "max_reuse_prop":max(float(z["max_control_reuse_prop"]) for z in cq if z["max_control_reuse_prop"]!=""),
  "minimum_inverse_simpson_effective_support":min(float(z["effective_support_inverse_simpson"]) for z in cq if z["effective_support_inverse_simpson"]!=""),
  "reuse_status":"DIAGNOSTIC"}
]
write_tsv(OUT/"MAR06_matching_gate.tsv",list(gate[0].keys()),gate)

sd=OUT/"source_data"; sd.mkdir(exist_ok=True)
for s,d in [
 ("MAR06_raw26_candidate_support.tsv","MAR06_SOURCE_RAW26_CANDIDATE_SUPPORT.tsv"),
 ("MAR06_raw26_reuse_summary.tsv","MAR06_SOURCE_RAW26_REUSE_SUMMARY.tsv"),
 ("MAR06_raw26_set_integrity.tsv","MAR06_SOURCE_RAW26_SET_INTEGRITY.tsv"),
 ("MAR06_raw26_balance.tsv","MAR06_SOURCE_RAW26_BALANCE.tsv"),
 ("MAR06_common_matching_qc.tsv","MAR06_SOURCE_COMMON_MATCHING_QC.tsv"),
 ("MAR06_common_control_reuse.tsv","MAR06_SOURCE_COMMON_CONTROL_REUSE.tsv"),
 ("MAR06_common_candidate_pool_qc.tsv","MAR06_SOURCE_COMMON_CANDIDATE_POOL_QC.tsv"),
 ("MAR06_matching_gate.tsv","MAR06_SOURCE_MATCHING_GATE.tsv"),
]:
    shutil.copy2(OUT/s,sd/d)
summary=[
 {"metric":"MAR06_status","value":"PASS","note":"Raw26 K50, K20 sensitivity, and common-side hard gates pass"},
 {"metric":"raw26_primary_B","value":"10000","note":"K=50 primary"},
 {"metric":"raw26_primary_pool_per_anchor","value":"50","note":"26 anchors"},
 {"metric":"raw26_sensitivity_B","value":"10000","note":"K=20 sensitivity"},
 {"metric":"raw26_sensitivity_pool_per_anchor","value":"20","note":"26 anchors"},
 {"metric":"common_targets","value":"20","note":"7 higher-order + 13 standalones"},
 {"metric":"common_B_per_target","value":"10000","note":""},
 {"metric":"reuse_concentration","value":"DIAGNOSTIC","note":"Not automatic failure criterion"},
 {"metric":"rare_constraint_mutability_sensitivity","value":"NOT_RUN_WITH_REASON","note":"Required frozen fields absent in WP01 model frame"},
 {"metric":"common_SNP_local_density_sensitivity","value":"NOT_RUN_WITH_REASON","note":"Local gene density absent in frozen gene table"},
 {"metric":"next_stage","value":"MAR07","note":""},
]
write_tsv(sd/"MAR06_SOURCE_WRITEUP_SUMMARY.tsv",list(summary[0].keys()),summary)

files=[p for p in OUT.rglob("*") if p.is_file() and p.name not in {"MAR06_checksums.sha256","MAR06_PASS.txt"}]
lock={"stage":"MAR06","status":"PASS","version":"1.0R2","timestamp_utc":utc(),
      "preanalysis_lock_sha256":sha256_file(OUT/"MAR06_preanalysis_lock.json"),
      "output_hashes":{str(p.relative_to(OUT)):sha256_file(p) for p in files}}
(OUT/"MAR06_analysis_lock.json").write_text(json.dumps(lock,indent=2),encoding="utf-8")
lines=[f"{sha256_file(p)}  {p.relative_to(MA)}" for p in sorted([p for p in OUT.rglob("*") if p.is_file() and p.name!="MAR06_checksums.sha256"])]
(OUT/"MAR06_checksums.sha256").write_text("\n".join(lines)+"\n",encoding="utf-8")
(OUT/"MAR06_PASS.txt").write_text(
 "status=PASS\nstage=MAR06_MATCHED_NULL_SUPPORT_REUSE_BALANCE_AUDIT\nversion=1.0R2\n"
 "raw26_primary_K50=PASS\nraw26_sensitivity_K20=PASS\ncommon_side=PASS\n"
 "reuse_concentration=DIAGNOSTIC\nnext_stage=MAR07\n",encoding="utf-8")
print("MAR06_POSTFLIGHT=PASS")
