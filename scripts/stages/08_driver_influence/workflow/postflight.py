from pathlib import Path
import csv, gzip, json, hashlib, os, math, statistics, subprocess, shutil, re, itertools, collections
from datetime import datetime, timezone

CH3=Path(os.environ.get("CH3_ROOT","/home/h3021/chapter3"))
MA=Path(os.environ.get("MA_ROOT","/home/h3021/chapter3/12_molecular_autism_revision"))
PKG=Path(__file__).resolve().parent.parent
OUT=MA/"08_driver_influence"/"MAR08_failure_localization_influence_v1p1"
MAR02_SD=MA/"02_raw26_domains"/"MAR02_raw26_compression_v1"/"source_data"
MAR03_SD=MA/"03_common_retest"/"MAR03_common_variant_retest_v1.0R2"/"source_data"
MAR07_SD=MA/"07_chain_calibration"/"MAR07_chain_calibration_v1.0R4p1"/"source_data"
MAR07_RUN=MA/"07_chain_calibration"/"MAR07_chain_calibration_v1.0R4p1"
PGC=CH3/"06_magma_pgc"/"MDV6_v1_1"
SPARK=CH3/"07_magma_spark"/"MDV7_v1_1"
CONTRACT=PKG/"config"/"MAR08_REVISED_ANALYSIS_CONTRACT_v1.0.json"
PYTHON_BIN=Path(os.environ.get("MAR08_PYTHON",str(MA/"envs"/"mar02_mdv5_legacy"/"bin"/"python")))

def utc(): return datetime.now(timezone.utc).isoformat()
def ensure_dirs():
    for d in ["state","logs","work","work/sets","work/magma","analysis","source_data"]:
        (OUT/d).mkdir(parents=True,exist_ok=True)
def sha256_file(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    return h.hexdigest()
def open_text(p):
    return gzip.open(p,"rt",encoding="utf-8-sig",errors="replace",newline="") if str(p).endswith(".gz") else open(p,"rt",encoding="utf-8-sig",errors="replace",newline="")
def read_tsv(p):
    with open_text(p) as f:
        r=csv.DictReader(f,delimiter="\t"); return r.fieldnames or [],list(r)
def write_tsv(p,fields,rows):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True)
    opener=gzip.open if str(p).endswith(".gz") else open
    with opener(p,"wt" if str(p).endswith(".gz") else "w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore"); w.writeheader(); w.writerows(rows)
def bh(pvals):
    n=len(pvals)
    if not n:return []
    order=sorted(range(n),key=lambda i:(float(pvals[i]),i)); q=[1.0]*n; cur=1.0
    for rev,i in enumerate(reversed(order),1):
        rank=n-rev+1; val=min(1.0,float(pvals[i])*n/rank); cur=min(cur,val); q[i]=cur
    return q
def qtile(v,p):
    v=sorted(float(x) for x in v); n=len(v)
    if not n:return float("nan")
    pos=(n-1)*p; lo=int(math.floor(pos)); hi=int(math.ceil(pos)); f=pos-lo
    return v[lo]*(1-f)+v[hi]*f
def corr(a,b):
    a=[float(x) for x in a]; b=[float(x) for x in b]
    ma=statistics.mean(a); mb=statistics.mean(b)
    den=(sum((x-ma)**2 for x in a)*sum((y-mb)**2 for y in b))**0.5
    return sum((x-ma)*(y-mb) for x,y in zip(a,b))/den if den else float("nan")
def split_genes(s): return [x.strip() for x in str(s or "").split(";") if x.strip()]
def sign(x):
    x=float(x); return 1 if x>0 else -1 if x<0 else 0
def parse_token(p):
    d={}
    for line in Path(p).read_text(encoding="utf-8",errors="replace").splitlines():
        if "=" in line:
            k,v=line.split("=",1); d[k]=v
    return d
def load_contract(): return json.loads(CONTRACT.read_text(encoding="utf-8"))
def derive_driver_sets(rep_rows,unit_summary):
    unit_drivers=collections.defaultdict(set)
    rep_driver_inc=collections.defaultdict(lambda:collections.Counter())
    standalone_drivers={}
    for z in rep_rows:
        genes=split_genes(z.get("raw26_drivers",""))
        if z.get("assignment_type")=="HIGHER_ORDER_REPRESENTATIVE":
            u=z["unit_id"]
            for g in genes:
                unit_drivers[u].add(g); rep_driver_inc[u][g]+=1
        elif z.get("assignment_type")=="STANDALONE_REPRESENTATIVE":
            sid=z["standalone_id"]
            standalone_drivers[sid]=set(genes)
            for g in genes: rep_driver_inc[sid][g]+=1
    # enforce frozen higher-order driver counts
    sm={z["unit_id"]:z for z in unit_summary}
    for u in sm:
        exp=int(sm[u]["raw26_driver_n"]); obs=len(unit_drivers[u])
        if exp!=obs: raise RuntimeError(f"driver count mismatch {u}: frozen={exp} derived={obs}")
    return dict(unit_drivers),standalone_drivers,rep_driver_inc
def build_scenarios(unit_drivers,unit_summary):
    sm={z["unit_id"]:z for z in unit_summary}
    allgenes=sorted(set().union(*[set(v) for v in unit_drivers.values()]))
    rows=[]
    for g in allgenes:
        origin=sorted(u for u,gs in unit_drivers.items() if g in gs)
        rows.append({"scenario_type":"LEAVE_ONE","scenario_id":f"L1_{g}","omitted_genes":g,"origin_units":";".join(origin),"eligibility_units":";".join(origin)})
    pair_units=collections.defaultdict(set)
    for u,gs in unit_drivers.items():
        n=int(sm[u]["raw26_driver_n"])
        if 3<=n<=10:
            for a,b in itertools.combinations(sorted(gs),2): pair_units[(a,b)].add(u)
    for (a,b),us in sorted(pair_units.items()):
        rows.append({"scenario_type":"LEAVE_TWO","scenario_id":f"L2_{a}__{b}","omitted_genes":f"{a};{b}","origin_units":";".join(sorted(us)),"eligibility_units":";".join(sorted(us))})
    return rows

ensure_dirs(); checks=[]
def ck(i,ok,o,e,d=""): checks.append({"check_id":i,"status":"PASS" if ok else "HOLD","observed":o,"expected":e,"detail":d})
resolved=json.loads((OUT/"MAR08_resolved_inputs.json").read_text()); lock=json.loads((OUT/"MAR08_preanalysis_lock.json").read_text())
req=[OUT/"analysis"/"MAR08_CHAIN_FAILURE_LOCALIZATION.tsv",OUT/"analysis"/"MAR08_OBSERVED_VS_NULL_BURDEN.tsv",OUT/"analysis"/"MAR08_DRIVER_CONCENTRATION.tsv",OUT/"analysis"/"MAR08_BASELINE_REPLAY.tsv",OUT/"analysis"/"MAR08_LEAVE_ONE_DRIVER.tsv",OUT/"analysis"/"MAR08_LEAVE_TWO_DRIVER.tsv.gz"]
for p in req: ck("exists_"+p.name,p.is_file(),str(p),"exists")
# recurrence conditional output
rec=OUT/"analysis"/"MAR08_NULL_PATHWAY_RECURRENCE.tsv"; nr=OUT/"analysis"/"MAR08_NULL_PATHWAY_RECURRENCE_NOT_RUN.txt"; ck("module_F_output",rec.is_file() or nr.is_file(),"recurrence" if rec.is_file() else "NOT_RUN" if nr.is_file() else "missing","one valid output")
# input immutability
for k,exp in lock["input_hashes"].items():
    p=Path(resolved[k]); act=sha256_file(p) if p.is_file() else "MISSING"; ck("input_immutable_"+k,act==exp,act,exp)
if all(p.is_file() for p in req):
    _,rp=read_tsv(OUT/"analysis"/"MAR08_BASELINE_REPLAY.tsv"); ck("baseline_replay_n",len(rp)==14,len(rp),14); ck("baseline_replay_pass",all(z["status"]=="PASS" for z in rp),sum(z["status"]=="PASS" for z in rp),14)
    _,l1=read_tsv(OUT/"analysis"/"MAR08_LEAVE_ONE_DRIVER.tsv"); _,l2=read_tsv(OUT/"analysis"/"MAR08_LEAVE_TWO_DRIVER.tsv.gz")
    ck("leave1_rows",len(l1)==154,len(l1),154); ck("leave2_rows",len(l2)==504,len(l2),504)
    k1={(z["dataset"],z["scenario_id"],z["target_id"]) for z in l1}; k2={(z["dataset"],z["scenario_id"],z["target_id"]) for z in l2}; ck("leave1_unique",len(k1)==len(l1),len(k1),len(l1)); ck("leave2_unique",len(k2)==len(l2),len(k2),len(l2))
write_tsv(OUT/"MAR08_postflight_checks.tsv",["check_id","status","observed","expected","detail"],checks)
if any(z["status"]!="PASS" for z in checks): raise SystemExit("HOLD: MAR08 postflight failures")
# influence summary
_,base=read_tsv(Path(resolved["MAR03_PRIMARY"])); bm={(z["dataset"],z["target_id"]):z for z in base}; _,l1=read_tsv(OUT/"analysis"/"MAR08_LEAVE_ONE_DRIVER.tsv"); _,l2=read_tsv(OUT/"analysis"/"MAR08_LEAVE_TWO_DRIVER.tsv.gz")
summary=[]
for ds in ["PGC","SPARK"]:
  for t in sorted({z["target_id"] for z in base}):
    a1=[z for z in l1 if z["dataset"]==ds and z["target_id"]==t and z["affected_by_omission"]=="YES"]; a2=[z for z in l2 if z["dataset"]==ds and z["target_id"]==t and z["affected_by_omission"]=="YES"]
    sig1=sum(z["scenario_primary_positive"]=="YES" for z in a1); sig2=sum(z["scenario_primary_positive"]=="YES" for z in a2)
    summary.append({"dataset":ds,"target_id":t,"baseline_beta":bm[(ds,t)]["beta"],"baseline_q_BH":bm[(ds,t)]["q_BH"],"leave1_affected_n":len(a1),"leave1_max_abs_delta_beta":max([abs(float(z["delta_beta_vs_frozen"])) for z in a1],default=0),"leave1_sign_flip_n":sum(z["sign_flip"]=="YES" for z in a1),"leave1_q_lt_005_n":sig1,"leave1_classification_switch_n":sum(z["classification_switch"]=="YES" for z in a1),"leave2_affected_n":len(a2),"leave2_max_abs_delta_beta":max([abs(float(z["delta_beta_vs_frozen"])) for z in a2],default=0),"leave2_sign_flip_n":sum(z["sign_flip"]=="YES" for z in a2),"leave2_q_lt_005_n":sig2,"leave2_classification_switch_n":sum(z["classification_switch"]=="YES" for z in a2),"revised_label":"INFLUENCE_SENSITIVE_SENSITIVITY_ONLY_NOT_RESCUE" if sig1 or sig2 else "COMMON_SIDE_NULL_ROBUST_TO_DRIVER_OMISSION"})
write_tsv(OUT/"MAR08_INFLUENCE_SUMMARY.tsv",list(summary[0].keys()),summary)
any_l1=any(int(z["leave1_q_lt_005_n"])>0 for z in summary); any_l2=any(int(z["leave2_q_lt_005_n"])>0 for z in summary)
rec_status="MDV3_RECURRENCE_RUN_TIER1_NOT_RUN" if rec.is_file() else "NOT_RUN_WITH_REASON_PATHWAY_LEVEL_NULL_STATE_NOT_FROZEN"
dec=[
 {"decision_domain":"MAR07_override","status":"FAIL_REMAINS_BINDING","manuscript_implication":"MAR08 cannot restore confirmatory Tier1 status."},
 {"decision_domain":"leave_one_driver_common_side","status":"INFLUENCE_SENSITIVE" if any_l1 else "NULL_ROBUST","manuscript_implication":"Any omission-created q<0.05 is sensitivity-only and cannot rescue MAR03/MAR07." if any_l1 else "Common-side null classification was robust to omission of individual frozen rare-anchor drivers."},
 {"decision_domain":"leave_two_driver_common_side","status":"INFLUENCE_SENSITIVE" if any_l2 else "NO_PRIMARY_RESCUE","manuscript_implication":"Leave-two results are limited diagnostics only."},
 {"decision_domain":"observed_vs_null_burden","status":"POST_HOC_OMNIBUS_EXPLORATORY","manuscript_implication":"Burden tails may contextualize the pathway layer but cannot validate individual pathways."},
 {"decision_domain":"null_pathway_recurrence","status":rec_status,"manuscript_implication":"Only pre-existing pathway-level sidecars were used; no global-null rerun was permitted."},
 {"decision_domain":"MAR09_release","status":"YES_REVISED_CLAIM_HIERARCHY","manuscript_implication":"Proceed to MAR09 only with Tier1 as exploratory/prioritization output, common-variant conclusion unchanged, and MAR07 failure explicit."}
]
write_tsv(OUT/"MAR08_MANUSCRIPT_DECISION.tsv",list(dec[0].keys()),dec)
# canonical source data
for name in ["MAR08_CHAIN_FAILURE_LOCALIZATION.tsv","MAR08_OBSERVED_VS_NULL_BURDEN.tsv","MAR08_DRIVER_CONCENTRATION.tsv","MAR08_LEAVE_ONE_DRIVER.tsv","MAR08_BASELINE_REPLAY.tsv"]:
    shutil.copy2(OUT/"analysis"/name,OUT/"source_data"/("MAR08_SOURCE_"+name.replace("MAR08_","")))
shutil.copy2(OUT/"analysis"/"MAR08_LEAVE_TWO_DRIVER.tsv.gz",OUT/"source_data"/"MAR08_SOURCE_LEAVE_TWO_DRIVER.tsv.gz")
if rec.is_file(): shutil.copy2(rec,OUT/"source_data"/"MAR08_SOURCE_NULL_PATHWAY_RECURRENCE.tsv")
else: shutil.copy2(nr,OUT/"source_data"/"MAR08_SOURCE_NULL_PATHWAY_RECURRENCE_NOT_RUN.txt")
shutil.copy2(OUT/"MAR08_INFLUENCE_SUMMARY.tsv",OUT/"source_data"/"MAR08_SOURCE_INFLUENCE_SUMMARY.tsv"); shutil.copy2(OUT/"MAR08_MANUSCRIPT_DECISION.tsv",OUT/"source_data"/"MAR08_SOURCE_MANUSCRIPT_DECISION.tsv")
# locks/checksums
files=[p for p in OUT.rglob("*") if p.is_file() and p.name not in {"MAR08_checksums.sha256","MAR08_PASS.txt"}]
(OUT/"MAR08_analysis_lock.json").write_text(json.dumps({"stage":"MAR08_FAILURE_LOCALIZATION_AND_INFLUENCE","version":"1.0R1","timestamp_utc":utc(),"technical_status":"PASS","MAR07_override":"FAIL_REMAINS","output_hashes":{str(p.relative_to(OUT)):sha256_file(p) for p in files}},indent=2),encoding="utf-8")
allfiles=[p for p in OUT.rglob("*") if p.is_file() and p.name!="MAR08_checksums.sha256"]
(OUT/"MAR08_checksums.sha256").write_text("\n".join(f"{sha256_file(p)}  {p.relative_to(MA)}" for p in sorted(allfiles))+"\n",encoding="utf-8")
(OUT/"MAR08_PASS.txt").write_text("status=PASS\nstage=MAR08_FAILURE_LOCALIZATION_AND_INFLUENCE\ntechnical_status=PASS\nscientific_role=EXPLORATORY_DIAGNOSTIC_STATISTICAL_REASSESSMENT\nMAR07_status=FAIL_REMAINS_BINDING\nMAR09_release=YES_REVISED_CLAIM_HIERARCHY\n",encoding="utf-8")
print("MAR08_POSTFLIGHT=PASS"); print("MAR09_RELEASE=YES_REVISED_CLAIM_HIERARCHY")
