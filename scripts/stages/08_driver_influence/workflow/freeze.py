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

import sys
ensure_dirs(); pre=OUT/"MAR08_preflight_report.tsv"
if not pre.is_file(): raise SystemExit("HOLD: preflight first")
_,pr=read_tsv(pre)
if any(z["status"]!="PASS" for z in pr): raise SystemExit("HOLD: preflight not PASS")
resolved=json.loads((OUT/"MAR08_resolved_inputs.json").read_text())
sys.path.insert(0,str(PKG/"vendor")); import mar03_common as m3
_,us=read_tsv(Path(resolved["MAR02_UNIT_SUMMARY"])); _,rr=read_tsv(Path(resolved["MAR02_REPRESENTATIVES"])); _,mem=read_tsv(Path(resolved["MAR02_UNIT_MEMBERSHIP"])); _,preview=read_tsv(OUT/"MAR08_SCENARIO_PREVIEW.tsv")
ud,sd,inc=derive_driver_sets(rr,us)
units=sorted(z["unit_id"] for z in us)
# frozen expanded symbol membership
members={u:set() for u in units}; hgnc={}
for z in mem:
    if z.get("Expanded_flag")=="1" and z["unit_id"] in members:
        g=(z.get("gene") or "").strip().upper(); hg=(z.get("HGNC_id") or "").strip().upper()
        if g: members[z["unit_id"]].add(g); hgnc[g]=hg
sm={z["unit_id"]:z for z in us}
for u in units:
    if len(members[u])!=int(sm[u]["Expanded_n"]): raise SystemExit(f"HOLD: expanded n mismatch {u}")
write_tsv(OUT/"MAR08_FROZEN_PRIMARY_MEMBERSHIP.tsv",["unit_id","gene","HGNC_id","Expanded_flag"],[{"unit_id":u,"gene":g,"HGNC_id":hgnc.get(g,""),"Expanded_flag":1} for u in units for g in sorted(members[u])])
# scenario manifest baseline + frozen scenarios
sc=[{"scenario_type":"BASELINE","scenario_id":"BASELINE","omitted_genes":"","origin_units":"","eligibility_units":""}]+preview
write_tsv(OUT/"MAR08_FROZEN_SCENARIO_MANIFEST.tsv",["scenario_type","scenario_id","omitted_genes","origin_units","eligibility_units"],sc)
# map exact symbols via frozen gene maps and write MAGMA set files
maps={"PGC":m3.detect_map(Path(resolved["PGC_MAP"])),"SPARK":m3.detect_map(Path(resolved["SPARK_MAP"]))}
_,base=read_tsv(Path(resolved["MAR03_PRIMARY"])); bmap={(z["dataset"],z["target_id"]):z for z in base}
coverage=[]
set_hashes={}
for ds in ["PGC","SPARK"]:
    ddir=OUT/"work"/"sets"/ds; ddir.mkdir(parents=True,exist_ok=True)
    mp=maps[ds]
    for s in sc:
        omitted=set(split_genes(s["omitted_genes"]))
        sf=ddir/f"{s['scenario_id']}.sets"
        with open(sf,"w",encoding="utf-8",newline="\n") as f:
            for u in units:
                syms=sorted(members[u]-omitted)
                gids=[]; missing=[]
                for g in syms:
                    gid=mp["by_symbol"].get(g)
                    if not gid and hgnc.get(g): gid=mp["by_hgnc"].get(hgnc[g])
                    if gid:gids.append(gid)
                    else: missing.append(g)
                gids=sorted(set(gids),key=lambda x:(not str(x).isdigit(),str(x)))
                f.write(u+" "+" ".join(gids)+"\n")
                removed=len(members[u]&omitted)
                coverage.append({"dataset":ds,"scenario_type":s["scenario_type"],"scenario_id":s["scenario_id"],"target_id":u,"frozen_n":len(members[u]),"post_omission_symbol_n":len(syms),"mapped_n":len(gids),"missing_n":len(missing),"missing_genes":";".join(missing),"removed_symbol_n":removed,"affected_by_omission":"YES" if removed else "NO"})
                if s["scenario_id"]=="BASELINE":
                    exp=int(bmap[(ds,u)]["mapped_n"])
                    if len(gids)!=exp: raise SystemExit(f"HOLD: baseline mapped_n mismatch {ds}/{u}: {len(gids)} vs {exp}")
        set_hashes[f"{ds}/{sf.name}"]=sha256_file(sf)
write_tsv(OUT/"MAR08_FROZEN_SCENARIO_COVERAGE.tsv",["dataset","scenario_type","scenario_id","target_id","frozen_n","post_omission_symbol_n","mapped_n","missing_n","missing_genes","removed_symbol_n","affected_by_omission"],coverage)
lock={
 "stage":"MAR08_FAILURE_LOCALIZATION_AND_INFLUENCE","version":"1.0R1","timestamp_utc":utc(),
 "MAR07_override":"FAIL_REMAINS","analysis_role":"EXPLORATORY_DIAGNOSTIC_STATISTICAL_REASSESSMENT",
 "scenario_counts":{"baseline":1,"leave_one":sum(z["scenario_type"]=="LEAVE_ONE" for z in sc),"leave_two":sum(z["scenario_type"]=="LEAVE_TWO" for z in sc)},
 "leave_one_drivers":sorted(set(g for z in sc if z["scenario_type"]=="LEAVE_ONE" for g in split_genes(z["omitted_genes"]))),
 "leave_two_eligible_units":[z["unit_id"] for z in us if 3<=int(z["raw26_driver_n"])<=10],
 "model":"Exact MAR03 primary marginal MAGMA competitive interface; --model direction=greater; same w10 gene results; seven-test BH per dataset/scenario",
 "baseline_replay_required":"YES","baseline_numeric_tolerance":1e-8,
 "module_F_rule":"MDV3 recurrence only if 1000/1000 pre-existing mdv3.tsv sidecars; Tier1 pathway recurrence NOT_RUN unless frozen pathway-level Tier1 state exists",
 "hard_prohibitions":load_contract()["hard_prohibitions"],
 "input_hashes":resolved["hashes"],"scenario_set_hashes":set_hashes,
 "association_scenario_results_read":"NO"
}
(OUT/"MAR08_preanalysis_lock.json").write_text(json.dumps(lock,indent=2),encoding="utf-8")
print("MAR08_FREEZE=PASS"); print("SCENARIO_COUNTS=BASELINE:1 L1:11 L2:36"); print("HOLD_FOR_MANUAL_GATE")
