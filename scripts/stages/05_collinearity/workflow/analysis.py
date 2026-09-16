
from pathlib import Path
import csv, json, hashlib, os, math, itertools, shutil, re
from datetime import datetime, timezone

MA = Path(os.environ.get("MA_ROOT","/home/h3021/chapter3/12_molecular_autism_revision"))
CH3 = Path(os.environ.get("CH3_ROOT","/home/h3021/chapter3"))
PKG = Path(__file__).resolve().parent.parent
OUT = MA / "05_magma_collinearity" / "MAR05_collinearity_v1"
MAR03_OUT = MA / "03_common_retest" / "MAR03_common_variant_retest_v1.0R2"
PGC = CH3 / "06_magma_pgc" / "MDV6_v1_1"
SPARK = CH3 / "07_magma_spark" / "MDV7_v1_1"
PY = Path(os.environ.get("MAR05_PYTHON", str(MA/"envs"/"mar02_mdv5_legacy"/"bin"/"python")))

def utc():
    return datetime.now(timezone.utc).isoformat()

def sha256_file(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()

def read_tsv(p):
    with open(p,encoding="utf-8-sig",errors="replace",newline="") as f:
        r=csv.DictReader(f,delimiter="\t")
        return r.fieldnames or [], list(r)

def write_tsv(p,fields,rows):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True)
    with open(p,"w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

def parse_pass(p):
    d={}
    for line in Path(p).read_text(encoding="utf-8").splitlines():
        if "=" in line:
            k,v=line.split("=",1); d[k]=v
    return d

def pick(header,names):
    m={x.strip().lower():x for x in header}
    for n in names:
        if n.lower() in m: return m[n.lower()]
    return None

def resolve_mar04_pass():
    candidates=[
        MA/"MAR04_FREEZE_RELEASE_v1.0R1"/"manifest"/"MAR04_FREEZE_PASS.txt",
        MA/"04_effect_boundary"/"MAR04_effect_boundary_v1"/"MAR04_PASS.txt",
    ]
    good=[p for p in candidates if p.is_file()]
    if not good:
        raise SystemExit("HOLD: MAR04 PASS not found")
    return good[0]

import numpy as np
from collections import defaultdict
OUT.mkdir(parents=True,exist_ok=True)
lockp=OUT/"MAR05_preanalysis_lock.json"; tokp=OUT/"MAR05_STOPA_RELEASE.txt"
if not lockp.is_file() or not tokp.is_file(): raise SystemExit("HOLD: freeze/release required")
tok=parse_pass(tokp)
if tok.get("preanalysis_lock_sha256")!=sha256_file(lockp): raise SystemExit("HOLD: release-lock mismatch")

_,mem=read_tsv(MAR03_OUT/"MAR03_frozen_membership.tsv")
labels={}
_,lr=read_tsv(PKG/"config"/"MAR05_UNIT_LABELS.tsv")
for z in lr: labels[z["unit_id"]]=z["biological_label"]
standlabels={}
_,sr=read_tsv(PKG/"config"/"MAR05_STANDALONE_LABELS.tsv")
for z in sr: standlabels[z["standalone_id"]]=z["pathway_name"]

def gene_universe(path):
    h,rows=read_tsv(path)
    gid=pick(h,["MAGMA_GENE","GENE","gene_id"])
    if not gid: raise SystemExit(f"HOLD: cannot identify gene ID in {path}: {h}")
    return [str(z[gid]).strip() for z in rows if str(z.get(gid,"")).strip()]

def memberships(ds,fam,univ):
    U=set(univ); out=defaultdict(set)
    for z in mem:
        if z["dataset"]==ds and z["family"]==fam:
            g=str(z["magma_gene_id"]).strip()
            if g in U: out[z["target_id"]].add(g)
    return dict(out)

def phi_stats(A,B,n):
    n11=len(A&B); n10=len(A-B); n01=len(B-A); n00=n-n11-n10-n01
    union=n11+n10+n01
    jac=n11/union if union else float("nan")
    den=((n11+n10)*(n01+n00)*(n11+n01)*(n10+n00))**0.5
    phi=(n11*n00-n10*n01)/den if den else float("nan")
    return n11,union,jac,phi

pair_rows=[]; vif_rows=[]; condition_rows=[]; design_corr_rows=[]; stand_pair_rows=[]; reconcile=[]
for ds,path,result_file in [
    ("PGC",PGC/"PGC2019_gene_results_w10.tsv",MAR03_OUT/"MAR03_pgc_higher_order.tsv"),
    ("SPARK",SPARK/"SPARK_EUR_gene_results_w10.tsv",MAR03_OUT/"MAR03_spark_higher_order.tsv"),
]:
    univ=gene_universe(path); idx={g:i for i,g in enumerate(univ)}
    A=memberships(ds,"A",univ); B=memberships(ds,"B",univ)
    units=sorted(A); stands=sorted(B)
    if len(units)!=7 or len(stands)!=13:
        raise SystemExit(f"HOLD: {ds} targets A/B={len(units)}/{len(stands)} expected 7/13")
    X=np.zeros((len(univ),len(units)),float)
    for j,u in enumerate(units):
        for g in A[u]: X[idx[g],j]=1.0

    # Reconcile actual NGENES with frozen MAR03 results.
    _,rr=read_tsv(result_file)
    ngen={z["target_id"]:int(float(z["NGENES"])) for z in rr}
    for u in units:
        reconcile.append({"dataset":ds,"target_id":u,"design_mapped_n":len(A[u]),
                          "MAR03_NGENES":ngen.get(u,""),"match":"YES" if ngen.get(u)==len(A[u]) else "NO"})

    # Primary pairwise overlap and design correlation.
    for i,j in itertools.combinations(range(len(units)),2):
        ui,uj=units[i],units[j]
        inter,union,jac,phi=phi_stats(A[ui],A[uj],len(univ))
        pair_rows.append({"dataset":ds,"unit_i":ui,"label_i":labels.get(ui,ui),"unit_j":uj,"label_j":labels.get(uj,uj),
                          "n_i":len(A[ui]),"n_j":len(A[uj]),"intersect_n":inter,"union_n":union,"Jaccard":jac,"phi":phi})
        design_corr_rows.append({"dataset":ds,"unit_i":ui,"unit_j":uj,"correlation":phi})

    # VIF.
    for j,u in enumerate(units):
        y=X[:,j]
        Z=np.column_stack([np.ones(len(univ)),np.delete(X,j,axis=1)])
        coef, *_ = np.linalg.lstsq(Z,y,rcond=None)
        resid=y-Z@coef
        sse=float(np.sum(resid**2)); sst=float(np.sum((y-y.mean())**2))
        r2=1-sse/sst
        vif=float("inf") if 1-r2<=1e-12 else 1/(1-r2)
        cls="HIGH" if vif>10 else "CAUTION" if vif>=5 else "GREEN"
        vif_rows.append({"dataset":ds,"unit_id":u,"biological_label":labels.get(u,u),
                         "mapped_n":len(A[u]),"R2_on_other_units":r2,"VIF":vif,"VIF_class":cls})

    # Condition number.
    Z=(X-X.mean(axis=0))/X.std(axis=0,ddof=0)
    s=np.linalg.svd(Z,full_matrices=False,compute_uv=False)
    k=float(s.max()/s.min()) if s.min()>1e-12 else float("inf")
    kcls="HIGH" if k>100 else "CAUTION" if k>=30 else "GREEN"
    condition_rows.append({"dataset":ds,"n_genes":len(univ),"n_units":len(units),"condition_number":k,
                           "condition_class":kcls,"min_singular_value":float(s.min()),"max_singular_value":float(s.max()),
                           "singular_values":";".join(str(float(x)) for x in s)})

    # Standalone overlap (secondary only).
    for i,j in itertools.combinations(range(len(stands)),2):
        si,sj=stands[i],stands[j]
        inter,union,jac,phi=phi_stats(B[si],B[sj],len(univ))
        stand_pair_rows.append({"dataset":ds,"standalone_i":si,"label_i":standlabels.get(si,si),
                                "standalone_j":sj,"label_j":standlabels.get(sj,sj),
                                "n_i":len(B[si]),"n_j":len(B[sj]),"intersect_n":inter,"union_n":union,
                                "Jaccard":jac,"phi":phi})

write_tsv(OUT/"MAR05_membership_reconciliation.tsv",list(reconcile[0].keys()),reconcile)
write_tsv(OUT/"MAR05_unit_overlap_long.tsv",list(pair_rows[0].keys()),pair_rows)
write_tsv(OUT/"MAR05_design_correlation_long.tsv",list(design_corr_rows[0].keys()),design_corr_rows)
write_tsv(OUT/"MAR05_vif.tsv",list(vif_rows[0].keys()),vif_rows)
write_tsv(OUT/"MAR05_condition_number.tsv",list(condition_rows[0].keys()),condition_rows)
write_tsv(OUT/"MAR05_standalone_overlap_long.tsv",list(stand_pair_rows[0].keys()),stand_pair_rows)

# Wide matrices for plotting/readability.
def matrix_file(ds,rows,value_key,name):
    us=sorted({z["unit_i"] for z in rows if z["dataset"]==ds}|{z["unit_j"] for z in rows if z["dataset"]==ds})
    lookup={(z["unit_i"],z["unit_j"]):z[value_key] for z in rows if z["dataset"]==ds}
    table=[]
    for u in us:
        r={"unit_id":u}
        for v in us:
            if u==v: val=1.0
            else: val=lookup.get((u,v),lookup.get((v,u),""))
            r[v]=val
        table.append(r)
    write_tsv(OUT/name.format(dataset=ds),["unit_id"]+us,table)

for ds in ["PGC","SPARK"]:
    matrix_file(ds,pair_rows,"Jaccard","MAR05_{dataset}_jaccard_matrix.tsv")
    matrix_file(ds,pair_rows,"phi","MAR05_{dataset}_phi_matrix.tsv")

# Native coefficient covariance availability.
coef_rows=[]
for ds in ["PGC","SPARK"]:
    base=MAR03_OUT/"magma"/ds
    candidates=[]
    if base.is_dir():
        for p in base.rglob("*"):
            if p.is_file() and any(tok in p.name.lower() for tok in ["covariance","coef_corr","coefficient_corr","correlation_matrix"]):
                candidates.append(str(p))
    coef_rows.append({"dataset":ds,"available":"YES_REVIEW_REQUIRED" if candidates else "NO",
                      "status":"REVIEW_REQUIRED" if candidates else "NOT_AVAILABLE_FROM_NATIVE_OUTPUT",
                      "candidate_files":";".join(candidates),
                      "substitute_audit":"DESIGN_CORRELATION_PHI"})
write_tsv(OUT/"MAR05_coefficient_correlation.tsv",list(coef_rows[0].keys()),coef_rows)

# Gate.
gate_rows=[]
overall_high=False
for ds in ["PGC","SPARK"]:
    vr=[z for z in vif_rows if z["dataset"]==ds]
    cr=next(z for z in condition_rows if z["dataset"]==ds)
    pr=[z for z in pair_rows if z["dataset"]==ds]
    maxv=max(float(z["VIF"]) for z in vr); maxj=max(float(z["Jaccard"]) for z in pr); k=float(cr["condition_number"])
    high=(maxv>10 or k>100); overall_high=overall_high or high
    gate_rows.append({"dataset":ds,"max_VIF":maxv,"condition_number":k,"max_pairwise_Jaccard":maxj,
                      "hard_high_collinearity":"YES" if high else "NO",
                      "marginal_model_status":"PRIMARY",
                      "conditional_reporting_gate":"SUPPLEMENT_ONLY" if high else "MAIN_TEXT_SUPPORTIVE_ELIGIBLE",
                      "reason":"VIF>10 or condition>100" if high else "No hard collinearity threshold exceeded"})
gate_rows.append({"dataset":"COMBINED","max_VIF":max(float(z["VIF"]) for z in vif_rows),
                  "condition_number":max(float(z["condition_number"]) for z in condition_rows),
                  "max_pairwise_Jaccard":max(float(z["Jaccard"]) for z in pair_rows),
                  "hard_high_collinearity":"YES" if overall_high else "NO",
                  "marginal_model_status":"PRIMARY",
                  "conditional_reporting_gate":"SUPPLEMENT_ONLY" if overall_high else "MAIN_TEXT_SUPPORTIVE_ELIGIBLE",
                  "reason":"At least one dataset exceeds hard threshold" if overall_high else "No PGC/SPARK hard threshold exceeded"})
write_tsv(OUT/"MAR05_conditional_model_gate.tsv",list(gate_rows[0].keys()),gate_rows)

print("MAR05_ANALYSIS=PASS")
