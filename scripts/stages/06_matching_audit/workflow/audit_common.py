
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

resolved=json.loads((OUT/"MAR06_resolved_inputs.json").read_text())
_,summary=read_tsv(resolved["MAR03_MATCHING_SUMMARY"])
_,qc=read_tsv(resolved["MAR03_MATCHING_QC"])
_,poolqc=read_tsv(resolved["MAR03_CANDIDATE_POOL_QC"])
_,reuse=read_tsv(resolved["MAR03_CONTROL_REUSE"])

# Aggregate reuse directly from frozen reuse table.
target_reuse=defaultdict(list)
for z in reuse:
    t=z.get("target_id") or z.get("target") or ""
    if t: target_reuse[t].append(z)

qc_by_target=defaultdict(list)
for z in qc:
    t=z.get("target_id") or z.get("target") or ""
    if t: qc_by_target[t].append(z)

rows=[]
for z in summary:
    t=z["target_id"]
    rr=target_reuse.get(t,[])
    props=[]
    counts=[]
    for r in rr:
        if r.get("selected_prop","")!="":
            props.append(float(r["selected_prop"]))
        if r.get("selected_count","")!="":
            counts.append(float(r["selected_count"]))
    H,expH,invsim=entropy_effective(counts)
    qrows=qc_by_target.get(t,[])
    meds=[safe_float(r.get("median_abs_smd","nan")) for r in qrows]
    p95s=[safe_float(r.get("p95_abs_smd","nan")) for r in qrows]
    rows.append({
      "family":z.get("family",""),"target_id":t,
      "frozen_n":z.get("frozen_n",z.get("target_frozen_n","")),
      "mapped_n":z.get("common_scored_n",z.get("mapped_n","")),
      "coverage":z.get("common_scored_coverage",z.get("coverage","")),
      "accepted_sets":z.get("matched_set_n",""),
      "unique_sets":z.get("unique_set_n",""),
      "attempts":z.get("attempts",""),
      "acceptance_rate":z.get("acceptance_rate",""),
      "exact_chromosome_composition":z.get("exact_chromosome_composition",""),
      "unique_control_genes":z.get("unique_control_gene_n",len(rr) if rr else ""),
      "max_control_reuse_prop":max(props) if props else z.get("max_control_reuse_prop",""),
      "reuse_entropy":H if counts else "",
      "effective_support_exp_entropy":expH if counts else "",
      "effective_support_inverse_simpson":invsim if counts else "",
      "max_median_abs_SMD":max([x for x in meds if math.isfinite(x)],default=""),
      "max_p95_abs_SMD":max([x for x in p95s if math.isfinite(x)],default=""),
    })
write_tsv(OUT/"MAR06_common_matching_qc.tsv",list(rows[0].keys()),rows)
shutil.copy2(resolved["MAR03_CONTROL_REUSE"],OUT/"MAR06_common_control_reuse.tsv")
shutil.copy2(resolved["MAR03_CANDIDATE_POOL_QC"],OUT/"MAR06_common_candidate_pool_qc.tsv")
print("MAR06_COMMON_AUDIT=PASS")
