
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

import numpy as np

resolved=json.loads((OUT/"MAR06_resolved_inputs.json").read_text())
authority=json.loads((PKG/"config"/"MAR06_RAW26_AUTHORITY.json").read_text())
all_support=[]
reuse_all=[]
integrity=[]
balance_all=[]

def load(role):
    rr=resolved["RAW26"][role]
    _,pool=read_tsv(rr["candidate_pool"])
    _,assign=read_tsv(rr["assignments"])
    _,sets=read_tsv(rr["sets"])
    _,balance=read_tsv(rr["balance_long"])
    _,bal_summary=read_tsv(rr["balance_summary"])
    return rr,pool,assign,sets,balance,bal_summary

def audit_role(role,label,K,B):
    rr,pool,assign,sets,balance,bal_summary=load(role)
    anchors=sorted({z["anchor_gene"] for z in pool})
    pool_by=defaultdict(list)
    for z in pool: pool_by[z["anchor_gene"]].append(z)
    ass_by_rep=defaultdict(list); sel=Counter(); sel_anchor=Counter()
    for z in assign:
        rep=str(z["replicate_id"]); a=z["anchor_gene"]; g=z["control_gene"]
        ass_by_rep[rep].append((a,g)); sel[g]+=1; sel_anchor[(a,g)]+=1
    set_by_rep={}
    set_hashes=[]
    for z in sets:
        rep=str(z["replicate_id"])
        genes=tuple(sorted([g for g in z["gene_set"].split(";") if g]))
        set_by_rep[rep]=genes
        set_hashes.append(z["set_sha256"])
    membership_mismatch=0; within_dup=0; assignment_anchor_dup=0
    for rep,pairs in ass_by_rep.items():
        controls=[g for a,g in pairs]
        if len(controls)!=len(set(controls)): within_dup+=1
        if len({a for a,g in pairs})!=26: assignment_anchor_dup+=1
        if tuple(sorted(controls))!=set_by_rep.get(rep,()): membership_mismatch+=1

    # Candidate support/reuse per anchor.
    for a in anchors:
        prow=pool_by[a]
        cands=sorted({z["candidate_gene"] for z in prow})
        counts=[sel_anchor[(a,c)] for c in cands]
        H,expH,invsim=entropy_effective(counts)
        reuse_all.append({
          "analysis":label,"anchor_gene":a,"n_candidates":len(cands),
          "selected_unique_candidates":sum(x>0 for x in counts),
          "max_reuse_prop":max(counts)/B if counts else "",
          "reuse_entropy":H,"effective_support_exp_entropy":expH,
          "effective_support_inverse_simpson":invsim,
          "accepted_sets":len(ass_by_rep),"unique_set_hashes":len(set(set_hashes)),
          "within_set_duplicate_sets":within_dup,"assignment_set_membership_mismatch":membership_mismatch
        })
        for z in prow:
            c=z["candidate_gene"]; cnt=sel_anchor[(a,c)]
            all_support.append({
              "analysis":label,"anchor_gene":a,"candidate_gene":c,
              "rank":z.get("rank",""),"distance":z.get("distance",""),
              "anchor_length_quintile":z.get("anchor_length_quintile",""),
              "candidate_length_quintile":z.get("candidate_length_quintile",""),
              "expansion_rule":f"{label}_{K}",
              "selected_count":cnt,"selected_prop":cnt/B
            })

    # Balance by covariate from long file.
    bycov=defaultdict(list)
    for z in balance:
        bycov[z["covariate"]].append(abs(float(z["SMD"])))
    for cov,vals in sorted(bycov.items()):
        arr=np.array(vals,float)
        med=float(np.median(arr)); p95=float(np.quantile(arr,0.95)); mx=float(np.max(arr))
        balance_all.append({
          "analysis":label,"covariate":cov,"n_replicates":len(vals),
          "median_abs_SMD":med,"p95_abs_SMD":p95,"max_abs_SMD":mx,
          "median_gate":"PASS" if med<=0.15 else "FAIL",
          "p95_gate":"PASS" if p95<=0.25 else "FAIL",
          "status":"PASS" if med<=0.15 and p95<=0.25 else "HOLD"
        })

    # Compare to frozen summary file if columns are available.
    summary_match="NOT_CHECKED"
    if bal_summary:
        smap={z["covariate"]:z for z in bal_summary if z.get("covariate")}
        diffs=[]
        for z in [x for x in balance_all if x["analysis"]==label]:
            s=smap.get(z["covariate"])
            if s:
                for col in ["median_abs_SMD","p95_abs_SMD","max_abs_SMD"]:
                    if s.get(col,"")!="":
                        diffs.append(abs(float(s[col])-float(z[col])))
        summary_match="PASS" if diffs and max(diffs)<=1e-12 else ("NOT_COMPARABLE" if not diffs else "FAIL")

    integ={
      "analysis":label,"K":K,"expected_B":B,"anchor_n":len(anchors),
      "candidate_pool_rows":len(pool),"candidate_pool_per_anchor_min":min(len(pool_by[a]) for a in anchors),
      "candidate_pool_per_anchor_max":max(len(pool_by[a]) for a in anchors),
      "assignment_rows":len(assign),"accepted_replicates":len(ass_by_rep),
      "set_rows":len(sets),"unique_set_hashes":len(set(set_hashes)),
      "within_set_duplicate_sets":within_dup,"assignment_anchor_duplicate_or_missing_sets":assignment_anchor_dup,
      "assignment_set_membership_mismatch":membership_mismatch,"balance_rows":len(balance),
      "balance_summary_replay":summary_match
    }
    integrity.append(integ)

audit_role("primary_K50","PRIMARY_K50",50,10000)
audit_role("sensitivity_K20","SENSITIVITY_K20",20,10000)

write_tsv(OUT/"MAR06_raw26_candidate_support.tsv",list(all_support[0].keys()),all_support)
write_tsv(OUT/"MAR06_raw26_reuse_summary.tsv",list(reuse_all[0].keys()),reuse_all)
write_tsv(OUT/"MAR06_raw26_set_integrity.tsv",list(integrity[0].keys()),integrity)
write_tsv(OUT/"MAR06_raw26_balance.tsv",list(balance_all[0].keys()),balance_all)

print("MAR06_RAW26_AUDIT=PASS")
for z in integrity:
    print(z["analysis"],"B=",z["accepted_replicates"],"unique=",z["unique_set_hashes"],
          "mismatch=",z["assignment_set_membership_mismatch"])
