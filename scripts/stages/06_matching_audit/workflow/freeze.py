
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

if not (OUT/"MAR06_preflight_report.tsv").is_file(): raise SystemExit("HOLD: preflight first")
_,pre=read_tsv(OUT/"MAR06_preflight_report.tsv")
if any(z["status"]!="PASS" for z in pre): raise SystemExit("HOLD: preflight not PASS")
resolved=json.loads((OUT/"MAR06_resolved_inputs.json").read_text())
authority=json.loads((PKG/"config"/"MAR06_RAW26_AUTHORITY.json").read_text())
lock={
 "stage":"MAR06_MATCHED_NULL_SUPPORT_REUSE_BALANCE_AUDIT","version":"1.0R2",
 "timestamp_utc":utc(),"audit_only_existing_frozen_matched_sets":"YES","regenerate_matched_sets":"NO",
 "raw26_primary":{"K":50,"B":10000,"candidate_pool_per_anchor":50,
   "artifacts":resolved["RAW26"]["primary_K50"],
   "balance_gate":authority["gates"]},
 "raw26_sensitivity":{"K":20,"B":10000,"candidate_pool_per_anchor":20,
   "artifacts":resolved["RAW26"]["sensitivity_K20"],
   "balance_gate":authority["gates"],
   "role":"PRESPECIFIED_SENSITIVITY"},
 "common_side":{"targets_expected":20,"B_expected":10000,"same_matched_sets_PGC_SPARK":"YES",
   "balance_gate":{"median_abs_SMD_max":0.15,"p95_abs_SMD_max":0.25}},
 "hard_fail":["accepted_sets_below_expected","duplicate_accepted_sets","assignment_set_membership_mismatch",
              "within_set_duplicate_genes","candidate_pool_size_mismatch","balance_gate_failure",
              "exact_chromosome_rule_violation"],
 "reuse_concentration_status":"DIAGNOSTIC_NOT_AUTOMATIC_FAILURE",
 "rare_side_constraint_mutability_sensitivity":"NOT_RUN_WITH_REASON_FIELDS_ABSENT_IN_FROZEN_WP01_MODEL_FRAME",
 "common_side_SNP_local_density_sensitivity":"NOT_RUN_WITH_REASON_LOCAL_GENE_DENSITY_ABSENT_IN_FROZEN_GENE_TABLE",
 "association_scores_used_for_matching_audit":"NO",
 "input_hashes":{}
}
for role in ["primary_K50","sensitivity_K20"]:
    for pth in resolved["RAW26"][role].values():
        p=Path(pth)
        if p.is_file(): lock["input_hashes"][str(p)]=sha256_file(p)
for key in ["MAR03_MATCHED_SETS","MAR03_MATCHING_SMD","MAR03_CANDIDATE_POOL_QC","MAR03_CONTROL_REUSE","MAR03_MATCHING_SUMMARY"]:
    p=Path(resolved[key])
    if p.is_file(): lock["input_hashes"][str(p)]=sha256_file(p)
p=OUT/"MAR06_preanalysis_lock.json"
p.write_text(json.dumps(lock,indent=2),encoding="utf-8")
print("MAR06_FREEZE=PASS")
print("RAW26_PRIMARY=K50_10000")
print("RAW26_SENSITIVITY=K20_10000")
print("HOLD_FOR_MANUAL_GATE")
