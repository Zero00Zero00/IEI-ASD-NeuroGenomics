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
SNAPSHOT_ROOT=PKG/"vendor"/"upstream_frozen"
SNAPSHOT_MANIFEST=PKG/"config"/"MAR08_UPSTREAM_FROZEN_SNAPSHOT_MANIFEST.json"

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

import sys, importlib.util
ensure_dirs(); checks=[]
def ck(i,ok,o,e="",d=""): checks.append({"check_id":i,"status":"PASS" if ok else "HOLD","observed":o,"expected":e,"detail":d})
contract=load_contract(); ck("contract_status",contract.get("status")=="FROZEN_BEFORE_ANALYSIS",contract.get("status"),"FROZEN_BEFORE_ANALYSIS")
exp_contract_sha="c0ed0c9b4d8ec94258666441c3766c67a0f932096eab1b65ed107f7d24e73578"; ck("contract_sha256",sha256_file(CONTRACT)==exp_contract_sha,sha256_file(CONTRACT),exp_contract_sha)
req={
 "MAR02_UNIT_SUMMARY":SNAPSHOT_ROOT/"MAR02"/"MAR02_SOURCE_UNIT_SUMMARY.tsv",
 "MAR02_UNIT_MEMBERSHIP":SNAPSHOT_ROOT/"MAR02"/"MAR02_SOURCE_UNIT_GENE_MEMBERSHIP.tsv",
 "MAR02_REPRESENTATIVES":SNAPSHOT_ROOT/"MAR02"/"MAR02_SOURCE_REPRESENTATIVES.tsv",
 "MAR02_STANDALONES":SNAPSHOT_ROOT/"MAR02"/"MAR02_SOURCE_STANDALONES.tsv",
 "MAR02_TIER1":SNAPSHOT_ROOT/"MAR02"/"MAR02_SOURCE_TIER1_PATHWAYS.tsv",
 "MAR03_PRIMARY":SNAPSHOT_ROOT/"MAR03"/"MAR03_SOURCE_PRIMARY_HIGHER_ORDER.tsv",
 "MAR07_NULL_DIST":MAR07_SD/"MAR07_SOURCE_FALSE_TIER1_DISTRIBUTION.tsv",
 "MAR07_GATE":MAR07_SD/"MAR07_SOURCE_CALIBRATION_GATE.tsv",
 "MAR04_SUMMARY":SNAPSHOT_ROOT/"MAR04"/"MAR04_SOURCE_WRITEUP_SUMMARY.tsv",
 "MAR05_SUMMARY":MA/"05_magma_collinearity"/"MAR05_collinearity_v1"/"source_data"/"MAR05_SOURCE_WRITEUP_SUMMARY.tsv",
 "MAR06_GATE":MA/"06_matched_null_audit"/"MAR06_matched_null_audit_v1.0R2"/"source_data"/"MAR06_SOURCE_MATCHING_GATE.tsv",
 "MAR02_CHECKSUMS":MA/"02_raw26_domains"/"MAR02_raw26_compression_v1"/"MAR02_checksums.sha256",
 "MAR03_CHECKSUMS":MA/"03_common_retest"/"MAR03_common_variant_retest_v1.0R2"/"MAR03_checksums.sha256",
 "MAR04_CHECKSUMS":MA/"04_effect_boundary"/"MAR04_effect_boundary_v1"/"MAR04_checksums.sha256",
 "MAR05_CHECKSUMS":MA/"05_magma_collinearity"/"MAR05_collinearity_v1"/"MAR05_checksums.sha256",
 "MAR06_CHECKSUMS":MA/"06_matched_null_audit"/"MAR06_matched_null_audit_v1.0R2"/"MAR06_checksums.sha256",
 "MAR07_CHECKSUMS":MAR07_RUN/"MAR07_checksums.sha256",
 "PGC_CHECKSUMS":PGC/"MDV6_checksums.sha256",
 "SPARK_CHECKSUMS":SPARK/"MDV7_checksums.sha256",
 "PGC_MAP":PGC/"MDV6_gene_id_map.tsv",
 "PGC_W10_RAW":PGC/"magma"/"PGC2019_w10.genes.raw",
 "PGC_PASS":PGC/"MDV6_PASS.txt",
 "SPARK_MAP":SPARK/"MDV7_gene_id_map.tsv",
 "SPARK_W10_RAW":SPARK/"magma"/"SPARK_EUR_w10.genes.raw",
 "SPARK_PASS":SPARK/"MDV7_PASS.txt",
 "MAR03_VENDOR_COMMON":PKG/"vendor"/"mar03_common.py",
 "MAR03_VENDOR_ANALYSIS":PKG/"vendor"/"mar03_analysis.py",
}
for k,p in req.items(): ck(k,p.is_file(),str(p),"exists")
# R1.1 authority transport snapshot checks.
snap=json.loads(SNAPSHOT_MANIFEST.read_text(encoding="utf-8"))
ck("snapshot_manifest_status",snap.get("status")=="FROZEN_TRANSPORT_SNAPSHOT",snap.get("status"),"FROZEN_TRANSPORT_SNAPSHOT")
contract_release_sha=contract.get("source_package_sha256",{})
for role in ["MAR02_UNIT_SUMMARY","MAR02_UNIT_MEMBERSHIP","MAR02_REPRESENTATIVES","MAR02_STANDALONES","MAR02_TIER1","MAR03_PRIMARY","MAR04_SUMMARY"]:
    p=req[role]
    if not p.is_file():
        continue
    name=p.name
    rec=snap["files"].get(name,{})
    exp=rec.get("sha256","")
    ck(f"{role}_snapshot_sha",bool(exp) and sha256_file(p)==exp,sha256_file(p),exp)
    release=rec.get("source_freeze_release","")
    release_sha=rec.get("source_freeze_release_sha256","")
    contract_sha=contract_release_sha.get(release,"")
    ck(f"{role}_release_provenance",bool(release_sha) and release_sha==contract_sha,release_sha,contract_sha,f"release={release}")
# frozen vendor hash
if req["MAR03_VENDOR_COMMON"].is_file(): ck("vendor_common_sha",sha256_file(req["MAR03_VENDOR_COMMON"])=="8f84aa9db69b98402ba1d72530785cb5c0384462f4e253c3bdea357bbeaa4844",sha256_file(req["MAR03_VENDOR_COMMON"]),"8f84aa9db69b98402ba1d72530785cb5c0384462f4e253c3bdea357bbeaa4844")
if req["MAR03_VENDOR_ANALYSIS"].is_file(): ck("vendor_analysis_sha",sha256_file(req["MAR03_VENDOR_ANALYSIS"])=="992f821d603bddacf4f879f6fbe4d6c532f0f1df733e223fdde6ed3265864570",sha256_file(req["MAR03_VENDOR_ANALYSIS"]),"992f821d603bddacf4f879f6fbe4d6c532f0f1df733e223fdde6ed3265864570")
# counts / scientific facts
if all(req[k].is_file() for k in ["MAR02_UNIT_SUMMARY","MAR02_STANDALONES","MAR02_TIER1","MAR02_REPRESENTATIVES"]):
    _,us=read_tsv(req["MAR02_UNIT_SUMMARY"]); _,sa=read_tsv(req["MAR02_STANDALONES"]); _,t1=read_tsv(req["MAR02_TIER1"]); _,rr=read_tsv(req["MAR02_REPRESENTATIVES"])
    ck("higher_order_n",len(us)==7,len(us),7); ck("standalone_n",len(sa)==13,len(sa),13); ck("tier1_n",len(t1)==39,len(t1),39); ck("representatives_n",len(rr)==30,len(rr),30)
    ud,sd,inc=derive_driver_sets(rr,us); sc=build_scenarios(ud,us)
    l1=[z for z in sc if z["scenario_type"]=="LEAVE_ONE"]; l2=[z for z in sc if z["scenario_type"]=="LEAVE_TWO"]
    ck("leave1_scenario_n",len(l1)==11,len(l1),11); ck("leave2_scenario_n",len(l2)==36,len(l2),36)
    elig=[z["unit_id"] for z in us if 3<=int(z["raw26_driver_n"])<=10]; ck("leave2_eligible_units",set(elig)=={"M04","M05"},";".join(sorted(elig)),"M04;M05")
    write_tsv(OUT/"MAR08_SCENARIO_PREVIEW.tsv",["scenario_type","scenario_id","omitted_genes","origin_units","eligibility_units"],sc)
if req["MAR03_PRIMARY"].is_file():
    _,p3=read_tsv(req["MAR03_PRIMARY"]); ck("MAR03_primary_rows",len(p3)==14,len(p3),14); ck("MAR03_primary_BH_positive",sum(float(z["q_BH"])<0.05 and float(z["beta"])>0 for z in p3)==0,sum(float(z["q_BH"])<0.05 and float(z["beta"])>0 for z in p3),0)
if req["MAR07_NULL_DIST"].is_file():
    _,nd=read_tsv(req["MAR07_NULL_DIST"]); ids={int(z["replicate_id"]) for z in nd}; rate=sum(int(z["n_Tier1"])>0 for z in nd)/len(nd) if nd else -1
    ck("MAR07_R",len(nd)==1000,len(nd),1000); ck("MAR07_ids",ids==set(range(1,1001)),len(ids),1000); ck("MAR07_rate",abs(rate-0.906)<1e-12,rate,0.906)
# conditional sidecar census
sidecars=sum((MAR07_RUN/"executor_work"/f"rep_{i:04d}"/"mdv3.tsv").is_file() for i in range(1,1001)); ck("MAR07_mdv3_sidecar_census",True,sidecars,"conditional only","1000 allows MDV3 recurrence; incomplete => NOT_RUN")
# resolve MAGMA using exact MAR03 vendor helper
sys.path.insert(0,str(PKG/"vendor")); import mar03_common as m3
try:
    magma,allm=m3.find_magma(); ck("MAGMA",True,str(magma),"executable")
except Exception as e:
    magma=None; ck("MAGMA",False,repr(e),"executable")
if PYTHON_BIN.is_file():
    cp=subprocess.run([str(PYTHON_BIN),"-c","import json,csv; print('OK')"],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True); ck("PYTHON_RUNTIME",cp.returncode==0,str(PYTHON_BIN),"validated",cp.stdout.strip())
else: ck("PYTHON_RUNTIME",False,str(PYTHON_BIN),"exists")

# verify frozen stage/source checksum manifests before using them
for label,target,manifest in [
 ("MAR05_summary_checksum",req["MAR05_SUMMARY"],req["MAR05_CHECKSUMS"]),
 ("MAR06_gate_checksum",req["MAR06_GATE"],req["MAR06_CHECKSUMS"]),
 ("MAR07_null_checksum",req["MAR07_NULL_DIST"],req["MAR07_CHECKSUMS"]),
 ("MAR07_gate_checksum",req["MAR07_GATE"],req["MAR07_CHECKSUMS"]),
 ("PGC_map_checksum",req["PGC_MAP"],req["PGC_CHECKSUMS"]),
 ("PGC_w10_checksum",req["PGC_W10_RAW"],req["PGC_CHECKSUMS"]),
 ("SPARK_map_checksum",req["SPARK_MAP"],req["SPARK_CHECKSUMS"]),
 ("SPARK_w10_checksum",req["SPARK_W10_RAW"],req["SPARK_CHECKSUMS"]),
]:
    if target.is_file() and manifest.is_file():
        st,exp=m3.verify_manifest_entry(manifest,target); ck(label,st=="PASS",st,"PASS",f"expected={exp}")
resolved={k:str(v) for k,v in req.items()}; resolved.update({"MAGMA_BIN":str(magma) if magma else "","PYTHON_BIN":str(PYTHON_BIN),"MAR07_MDV3_SIDECAR_N":sidecars,"UPSTREAM_FROZEN_SNAPSHOT_MANIFEST":str(SNAPSHOT_MANIFEST)})
resolved["hashes"]={k:sha256_file(v) for k,v in req.items() if v.is_file()}
(OUT/"MAR08_resolved_inputs.json").write_text(json.dumps(resolved,indent=2),encoding="utf-8")
write_tsv(OUT/"MAR08_preflight_report.tsv",["check_id","status","observed","expected","detail"],checks)
bad=[z for z in checks if z["status"]!="PASS"]
if bad: raise SystemExit(f"HOLD: MAR08 preflight failed ({len(bad)} checks); see MAR08_preflight_report.tsv")
print("MAR08_PREFLIGHT=PASS"); print(f"SCENARIOS=L1:{len(l1)} L2:{len(l2)}"); print(f"MAR07_MDV3_SIDECARS={sidecars}/1000")
