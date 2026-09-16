
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

import sys
ensure_dirs()
mode=sys.argv[1] if len(sys.argv)>1 else "pilot"
executor=PKG/"executor"/"MAR07_fullchain_executor.sh"
if not (OUT/"state"/"MAR07_FULLCHAIN_EXECUTOR_CERTIFIED.txt").is_file():
    raise SystemExit("HOLD: certification required")
lock=json.loads((OUT/"MAR07_preanalysis_lock.json").read_text())
cfg=lock["global_null"]; seed_base=int(cfg["seed_base"])

def run_rep(rep,seed):
    p=OUT/"calibration"/f"rep_{rep:04d}.json"
    if p.is_file(): return json.loads(p.read_text())
    cp=subprocess.run([str(executor),"--replicate-id",str(rep),"--seed",str(seed),"--out",str(p)],
                      stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    (OUT/"logs"/f"rep_{rep:04d}.log").write_text(cp.stdout,encoding="utf-8")
    if cp.returncode!=0 or not p.is_file(): raise SystemExit(f"HOLD: replicate {rep} failed")
    return json.loads(p.read_text())

if mode=="pilot": target=50
elif mode=="calibrate": target=250
else: raise SystemExit("pilot|calibrate")
res=[run_rep(r,seed_base+r) for r in range(1,target+1)]
if mode=="pilot":
    rows=[{"replicate_id":z["replicate_id"],"seed":z["seed"],"n_MDV2":z["n_MDV2"],"n_MDV3":z["n_MDV3"],"n_Tier1":z["n_Tier1"],
           "reference_B_MDV4":z.get("reference_B_MDV4",""),"runtime_sec":z["runtime_sec"]} for z in res]
    write_tsv(OUT/"calibration"/"MAR07_pilot_results.tsv",list(rows[0].keys()),rows)
    (OUT/"state"/"MAR07_PILOT_PASS.txt").write_text(
        f"status=PASS\nR=50\nmean_runtime_sec={sum(float(z['runtime_sec']) for z in res)/50}\n",encoding="utf-8")
    print("MAR07_PILOT=PASS"); raise SystemExit(0)

def summarize(rr):
    n=len(rr); k=sum(int(z["n_Tier1"])>0 for z in rr); rate=k/n; lo,hi=wilson_ci(k,n); half=(hi-lo)/2
    vals=[int(z["n_Tier1"]) for z in rr]; mean=sum(vals)/n; sd=statistics.stdev(vals) if n>1 else 0; se=sd/math.sqrt(n)
    return {"R":n,"any_false_tier1_n":k,"any_false_tier1_rate":rate,"MC_CI_low":lo,"MC_CI_high":hi,"MC_CI_halfwidth":half,
            "mean_false_tier1_n":mean,"mean_false_tier1_MC95_low":max(0,mean-1.95996398454*se),"mean_false_tier1_MC95_high":mean+1.95996398454*se}
summ=summarize(res)
for nxt in [500,1000]:
    if summ["MC_CI_halfwidth"]<=0.02: break
    for r in range(len(res)+1,nxt+1): res.append(run_rep(r,seed_base+r))
    summ=summarize(res)
rows=[{"replicate_id":z["replicate_id"],"seed":z["seed"],"n_MDV2":z["n_MDV2"],"n_MDV3":z["n_MDV3"],"n_Tier1":z["n_Tier1"],
       "reference_B_MDV4":z.get("reference_B_MDV4",""),"runtime_sec":z["runtime_sec"]} for z in res]
write_tsv(OUT/"calibration"/"MAR07_false_tier1_distribution.tsv",list(rows[0].keys()),rows)
write_tsv(OUT/"calibration"/"MAR07_chain_calibration_metrics.tsv",list(summ.keys()),[summ])
print("MAR07_CALIBRATION=PASS")
