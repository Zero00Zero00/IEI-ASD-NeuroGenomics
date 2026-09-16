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

ensure_dirs(); tok=OUT/"MA_MAR08_STOPA_RELEASE.txt"; lockp=OUT/"MAR08_preanalysis_lock.json"
if not tok.is_file() or not lockp.is_file(): raise SystemExit("HOLD: release required")
if parse_token(tok).get("preanalysis_lock_sha256")!=sha256_file(lockp): raise SystemExit("HOLD: release-lock mismatch")
resolved=json.loads((OUT/"MAR08_resolved_inputs.json").read_text())
_,dist=read_tsv(Path(resolved["MAR07_NULL_DIST"]));
if len(dist)!=1000: raise SystemExit("HOLD: MAR07 distribution is not R=1000")
# A: localization
loc=[]
vars=["n_MDV2","n_MDV3","n_Tier1"]
for v in vars:
    a=[int(z[v]) for z in dist]
    loc.append({"analysis":"COUNT_SUMMARY","variable_x":v,"variable_y":"","stratum":"ALL","n":len(a),"metric":"P_gt0","value":sum(x>0 for x in a)/len(a)})
    for metric,val in [("mean",statistics.mean(a)),("SD",statistics.stdev(a)),("median",statistics.median(a)),("Q1",qtile(a,.25)),("Q3",qtile(a,.75)),("P05",qtile(a,.05)),("P95",qtile(a,.95)),("min",min(a)),("max",max(a))]:
        loc.append({"analysis":"COUNT_SUMMARY","variable_x":v,"variable_y":"","stratum":"ALL","n":len(a),"metric":metric,"value":val})
for a,b in [("n_MDV2","n_MDV3"),("n_MDV2","n_Tier1"),("n_MDV3","n_Tier1"),("reference_B_MDV4","n_Tier1")]:
    loc.append({"analysis":"CORRELATION","variable_x":a,"variable_y":b,"stratum":"ALL","n":1000,"metric":"Pearson_r","value":corr([z[a] for z in dist],[z[b] for z in dist])})
for label,pred in [("B_LT_8000",lambda z:int(z["reference_B_MDV4"])<8000),("B_GE_8000",lambda z:int(z["reference_B_MDV4"])>=8000)]:
    sub=[z for z in dist if pred(z)]; a=[int(z["n_Tier1"]) for z in sub]
    for metric,val in [("P_gt0",sum(x>0 for x in a)/len(a)),("mean",statistics.mean(a)),("median",statistics.median(a)),("P05",qtile(a,.05)),("P95",qtile(a,.95)),("min",min(a)),("max",max(a))]:
        loc.append({"analysis":"B_STRATIFIED_TIER1","variable_x":"n_Tier1","variable_y":"reference_B_MDV4","stratum":label,"n":len(a),"metric":metric,"value":val})
write_tsv(OUT/"analysis"/"MAR08_CHAIN_FAILURE_LOCALIZATION.tsv",["analysis","variable_x","variable_y","stratum","n","metric","value"],loc)
# B: observed vs null burden
obs={"MDV2_primary_pass":("n_MDV2",18),"MDV3_BH_positive":("n_MDV3",205),"Tier1_BH":("n_Tier1",39)}; burden=[]
for stage,(col,o) in obs.items():
    a=[int(z[col]) for z in dist]; ge=sum(x>=o for x in a); sd=statistics.stdev(a); zpos=(o-statistics.mean(a))/sd if sd else float("nan")
    burden.append({"stage":stage,"null_column":col,"observed_count":o,"R":1000,"null_mean":statistics.mean(a),"null_SD":sd,"null_median":statistics.median(a),"null_P05":qtile(a,.05),"null_P95":qtile(a,.95),"null_ge_observed_n":ge,"P_emp_upper_plus1":(1+ge)/1001,"standardized_position_z":zpos,"interpretation":"POST_HOC_OMNIBUS_EXPLORATORY_NOT_INDIVIDUAL_PATHWAY_VALIDATION"})
write_tsv(OUT/"analysis"/"MAR08_OBSERVED_VS_NULL_BURDEN.tsv",list(burden[0].keys()),burden)
# C: driver concentration
_,us=read_tsv(Path(resolved["MAR02_UNIT_SUMMARY"])); _,rr=read_tsv(Path(resolved["MAR02_REPRESENTATIVES"])); _,sa=read_tsv(Path(resolved["MAR02_STANDALONES"]))
ud,sd,inc=derive_driver_sets(rr,us); sm={z["unit_id"]:z for z in us}; sam={z["standalone_id"]:z for z in sa}
rows=[]
for target,fam in [(u,"A_HIGHER_ORDER") for u in sorted(ud)]+[(s,"B_STANDALONE") for s in sorted(sd)]:
    ctr=inc[target]; total=sum(ctr.values()); drivers=sorted(ctr); weights={g:ctr[g]/total for g in drivers}; mx=max(ctr.values()) if ctr else 0; dom=sorted(g for g in drivers if ctr[g]==mx)
    neff=1/sum(p*p for p in weights.values()) if weights else float("nan"); cmax=max(weights.values()) if weights else float("nan")
    frozen_n=int(sm[target]["raw26_driver_n"]) if fam=="A_HIGHER_ORDER" else int(sam[target]["raw26_driver_n"])
    if frozen_n!=len(drivers): raise SystemExit(f"HOLD: driver_n mismatch {target}")
    rows.append({"target_id":target,"target_family":fam,"biological_label":sm[target]["biological_label"] if fam=="A_HIGHER_ORDER" else sam[target]["biological_label"],"representative_n":int(sm[target]["representative_count"]) if fam=="A_HIGHER_ORDER" else 1,"raw26_driver_n":len(drivers),"driver_genes":";".join(drivers),"dominant_driver_gene":";".join(dom),"dominant_tie_n":len(dom),"locus_dominated_flag":sm[target]["locus_dominated_flag"] if fam=="A_HIGHER_ORDER" else "NOT_FROZEN_FOR_STANDALONE","representative_pathway_driver_incidence":";".join(f"{g}:{ctr[g]}" for g in drivers),"driver_weights":";".join(f"{g}:{weights[g]:.6f}" for g in drivers),"N_eff":neff,"C_max":cmax,"interpretation":"DESCRIPTIVE_ONLY"})
write_tsv(OUT/"analysis"/"MAR08_DRIVER_CONCENTRATION.tsv",list(rows[0].keys()),rows)
print("MAR08_MODULES_A_B_C=PASS")
