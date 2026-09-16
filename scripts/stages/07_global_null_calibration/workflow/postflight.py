
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

ensure_dirs()
checks=[]
def ck(i,ok,o,e,d=""): checks.append({"check_id":i,"status":"PASS" if ok else "HOLD","observed":o,"expected":e,"detail":d})
req=[
 OUT/"state"/"MAR07_FULLCHAIN_EXECUTOR_CERTIFIED.txt",
 OUT/"observed"/"MAR07_observed_tier1_stability.tsv",
 OUT/"calibration"/"MAR07_chain_calibration_metrics.tsv",
 OUT/"calibration"/"MAR07_false_tier1_distribution.tsv"]
for p in req: ck("exists_"+p.name,p.is_file(),str(p),"exists")
if all(p.is_file() for p in req):
    _,tiers=read_tsv(req[1]); _,met=read_tsv(req[2]); _,dist=read_tsv(req[3])
    ck("Tier1_BH",sum(z["Tier1_BH"]=="YES" for z in tiers)==39,sum(z["Tier1_BH"]=="YES" for z in tiers),39)
    ck("Tier1_BY",sum(z["Tier1_BY"]=="YES" for z in tiers)==0,sum(z["Tier1_BY"]=="YES" for z in tiers),0)
    m=met[0]; R=int(float(m["R"])); ck("R",R in {250,500,1000},R,"250/500/1000"); ck("dist_R",len(dist)==R,len(dist),R)
write_tsv(OUT/"MAR07_postflight_checks.tsv",["check_id","status","observed","expected","detail"],checks)
bad=[z for z in checks if z["status"]!="PASS"]
if bad: raise SystemExit(f"HOLD: postflight {len(bad)} failures")
_,met=read_tsv(OUT/"calibration"/"MAR07_chain_calibration_metrics.tsv"); m=met[0]
rate=float(m["any_false_tier1_rate"]); hi=float(m["MC_CI_high"])
cal="FAIL" if hi>0.10 else "PASS" if rate<=0.05 else "CAUTION"
overall="FAIL" if cal=="FAIL" else "CAUTION"  # BY loses 39/39, so at least caution.
gate=[{"calibration_gate":cal,"any_false_tier1_rate":rate,"MC_CI_low":m["MC_CI_low"],"MC_CI_high":hi,
       "Tier1_BH_n":39,"Tier1_BY_n":0,"Tier1_lost_under_BY":39,"BY_stability":"SENSITIVE_REVIEW_REQUIRED",
       "overall_gate":overall,
       "manuscript_implication":"PAUSE_AND_REASSESS" if overall=="FAIL" else "REPORT_DEPENDENCY_SENSITIVITY_AND_EMPIRICAL_CALIBRATION"}]
write_tsv(OUT/"MAR07_CALIBRATION_GATE.tsv",list(gate[0].keys()),gate)

# Source data
for src,dst in [
 (OUT/"observed"/"MAR07_bh_by_sensitivity.tsv",OUT/"source_data"/"MAR07_SOURCE_BH_BY_SENSITIVITY.tsv"),
 (OUT/"observed"/"MAR07_observed_tier1_stability.tsv",OUT/"source_data"/"MAR07_SOURCE_TIER1_STABILITY.tsv"),
 (OUT/"observed"/"MAR07_observed_dependency_summary.tsv",OUT/"source_data"/"MAR07_SOURCE_OBSERVED_DEPENDENCY_SUMMARY.tsv"),
 (OUT/"calibration"/"MAR07_chain_calibration_metrics.tsv",OUT/"source_data"/"MAR07_SOURCE_CHAIN_CALIBRATION_METRICS.tsv"),
 (OUT/"calibration"/"MAR07_false_tier1_distribution.tsv",OUT/"source_data"/"MAR07_SOURCE_FALSE_TIER1_DISTRIBUTION.tsv"),
 (OUT/"MAR07_CALIBRATION_GATE.tsv",OUT/"source_data"/"MAR07_SOURCE_CALIBRATION_GATE.tsv"),
 (OUT/"certification"/"MAR07_FULLCHAIN_CERTIFICATION.tsv",OUT/"source_data"/"MAR07_SOURCE_FULLCHAIN_CERTIFICATION.tsv")]:
    shutil.copy2(src,dst)
summary=[
 {"metric":"MAR07_status","value":overall,"note":"BY sensitivity guarantees at least CAUTION unless policy amended prospectively"},
 {"metric":"Tier1_BH_n","value":"39","note":""},{"metric":"Tier1_BY_n","value":"0","note":""},
 {"metric":"global_null_R","value":m["R"],"note":""},{"metric":"P_any_false_Tier1","value":m["any_false_tier1_rate"],"note":""},
 {"metric":"MC_CI_high","value":m["MC_CI_high"],"note":"HOLD if >0.10"},
 {"metric":"next_stage","value":"MAR08" if overall!="FAIL" else "STATISTICAL_REASSESSMENT","note":""}]
write_tsv(OUT/"source_data"/"MAR07_SOURCE_WRITEUP_SUMMARY.tsv",list(summary[0].keys()),summary)

files=[p for p in OUT.rglob("*") if p.is_file() and p.name not in {"MAR07_checksums.sha256","MAR07_PASS.txt"}]
(OUT/"MAR07_analysis_lock.json").write_text(json.dumps({"stage":"MAR07","version":"1.0R4.1","status":overall,"timestamp_utc":utc(),
 "preanalysis_lock_sha256":sha256_file(OUT/"MAR07_preanalysis_lock.json"),
 "output_hashes":{str(p.relative_to(OUT)):sha256_file(p) for p in files}},indent=2),encoding="utf-8")
(OUT/"MAR07_checksums.sha256").write_text(
 "\n".join(f"{sha256_file(p)}  {p.relative_to(MA)}" for p in sorted([p for p in OUT.rglob("*") if p.is_file() and p.name!="MAR07_checksums.sha256"]))+"\n",encoding="utf-8")
(OUT/"MAR07_PASS.txt").write_text(
 f"status={overall}\nstage=MAR07_TARGETED_EMPIRICAL_CHAIN_CALIBRATION\nversion=1.0R4\ncalibration_gate={cal}\n"
 "BY_stability=SENSITIVE_REVIEW_REQUIRED\n"
 f"next_stage={'MAR08' if overall!='FAIL' else 'STATISTICAL_REASSESSMENT'}\n",encoding="utf-8")
print("MAR07_POSTFLIGHT="+overall)
