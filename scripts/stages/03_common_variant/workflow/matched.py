
from common import *
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from collections import defaultdict, Counter

ensure_dirs()
lockp=OUT/"MAR03_preanalysis_lock.json"; tokenp=OUT/"MAR03_STOPA_RELEASE.txt"
if not lockp.is_file() or not tokenp.is_file(): raise SystemExit("HOLD: freeze/release required")
lock=json.loads(lockp.read_text())
tok={}
for line in tokenp.read_text().splitlines():
    if "=" in line:
        k,v=line.split("=",1); tok[k]=v
if tok.get("preanalysis_lock_sha256")!=sha256_file(lockp):
    raise SystemExit("HOLD: release token mismatch")

mc=lock["matched_continuous"]
if mc["matching_method"]!="chromosome_stratified_exponential_tilt_rerandomization":
    raise SystemExit("HOLD: matching method not frozen to historical MDV8 Match-A1")

resolved=json.loads((OUT/"MAR03_resolved_inputs.json").read_text())
d=pd.read_csv(resolved["MDV8_common"],sep="\t")
d["MAGMA_GENE"]=d["MAGMA_GENE"].astype(str)
d=d.set_index("MAGMA_GENE",drop=False)

# Frozen target sets, identical PGC/SPARK IDs guaranteed at freeze.
_,mem=read_tsv(OUT/"MAR03_frozen_membership.tsv")
sets={}
for z in mem:
    if z["dataset"]=="PGC" and z["family"] in ("A","B"):
        sets.setdefault((z["family"],z["target_id"]),set()).add(str(z["magma_gene_id"]))
targets=sorted(sets, key=lambda x:(x[0],x[1]))

# Exclude all tested A+B hypothesis genes globally from control pool.
all_hyp=set().union(*[sets[k] for k in targets])
control=d[~d.MAGMA_GENE.isin(all_hyp)].copy()

feature_cols=mc["distance_variables"]
B=int(mc["B"]); seed=int(mc["seed"])
caliper=float(mc["per_set_abs_SMD_caliper"])
med_gate=float(mc["median_abs_SMD_max"])
p95_gate=float(mc["p95_abs_SMD_max"])
resid_gate=float(mc["tilt_mean_residual_max"])
max_factor=int(mc["max_set_generation_attempt_factor"])

def abs_smd_from_moments(target_mean,target_var,control_x):
    cm=np.mean(control_x,axis=0); cv=np.var(control_x,axis=0,ddof=1)
    den=np.sqrt((target_var+cv)/2.0); diff=np.abs(target_mean-cm)
    out=np.empty_like(diff)
    for j in range(len(diff)):
        if not np.isfinite(den[j]) or den[j]==0:
            out[j]=0.0 if np.isclose(diff[j],0.0) else np.inf
        else: out[j]=diff[j]/den[j]
    return out

def stable_weights(x,lam):
    z=np.asarray(x,float)@np.asarray(lam,float); z=z-np.max(z)
    w=np.exp(np.clip(z,-50.0,0.0)); s=float(w.sum())
    if not np.isfinite(s) or s<=0: raise RuntimeError("Invalid exponential-tilt weights")
    return w/s

def fit_tilt(target_mean,chr_counts,pools):
    total=int(sum(chr_counts.values()))
    def predicted(lam):
        acc=np.zeros(len(feature_cols),float)
        for ch,n in chr_counts.items():
            X=pools[ch]["X"]; w=stable_weights(X,lam)
            acc += int(n)*np.sum(w[:,None]*X,axis=0)
        return acc/total
    fit=least_squares(lambda lam: predicted(lam)-target_mean,
                      np.zeros(len(feature_cols),float),
                      max_nfev=500,xtol=1e-11,ftol=1e-11,gtol=1e-11)
    pred=predicted(fit.x); resid=pred-target_mean
    if not fit.success: raise RuntimeError(f"tilt solver failed: {fit.message}")
    if float(np.max(np.abs(resid)))>resid_gate:
        raise RuntimeError(f"tilt residual {np.max(np.abs(resid)):.6g} > {resid_gate}")
    return fit.x,pred,resid

matched_path=OUT/"matched"/"MAR03_matched_sets.tsv.gz"
smd_path=OUT/"matched"/"MAR03_matching_smd.tsv.gz"
pool_rows=[]; summary_rows=[]; qc_rows=[]; tilt_rows=[]; reuse=Counter()

with gzip.open(matched_path,"wt",encoding="utf-8") as sf, gzip.open(smd_path,"wt",encoding="utf-8") as mf:
    sf.write("family\ttarget_id\treplicate\tgenes\n")
    mf.write("family\ttarget_id\treplicate\tabs_smd_log_gene_length\tabs_smd_GC\tabs_smd_annotation\n")
    for ui,(fam,target) in enumerate(targets):
        target_ids=sorted(sets[(fam,target)] & set(d.index))
        frozen_n=len(sets[(fam,target)])
        coverage=len(target_ids)/frozen_n if frozen_n else 0
        if len(target_ids)<2: raise SystemExit(f"HOLD: {fam}/{target} <2 common-scored target genes")
        if coverage<0.8: raise SystemExit(f"HOLD: {fam}/{target} common-scored coverage {coverage:.3f}<0.8")
        t=d.loc[target_ids]
        tx=t[feature_cols].to_numpy(float)
        tmean=tx.mean(axis=0); tvar=tx.var(axis=0,ddof=1)
        chr_counts={int(ch):int(n) for ch,n in t.groupby("chr_mdv8").size().items()}
        pools={}
        for ch,n in chr_counts.items():
            cp=control[control.chr_mdv8==ch].copy()
            if len(cp)<n: raise SystemExit(f"HOLD: {target} chr {ch}: {len(cp)} controls < {n}")
            pools[ch]={"ids":cp.MAGMA_GENE.astype(str).to_numpy(),
                       "X":cp[feature_cols].to_numpy(float)}
            pool_rows.append({"family":fam,"target_id":target,"chromosome":ch,"target_chr_n":n,
                              "eligible_control_chr_n":len(cp),"global_hypothesis_union_excluded_n":len(all_hyp)})
        lam,pred,resid=fit_tilt(tmean,chr_counts,pools)
        for j,nm in enumerate(["log_gene_length","GC","annotation_degree"]):
            tilt_rows.append({"family":fam,"target_id":target,"covariate":nm,"lambda":float(lam[j]),
                              "target_mean":float(tmean[j]),"tilted_expected_mean":float(pred[j]),
                              "mean_residual":float(resid[j]),"abs_mean_residual":float(abs(resid[j]))})
        weighted={ch:stable_weights(p["X"],lam) for ch,p in pools.items()}
        rng=np.random.default_rng(seed+1009*ui)
        seen=set(); vals=[]; attempts=0; balrej=0; duprej=0; maxatt=B*max_factor
        while len(seen)<B and attempts<maxatt:
            attempts+=1; chosen=[]; chosen_x=[]
            for ch,n in chr_counts.items():
                p=pools[ch]; w=weighted[ch]
                idx=rng.choice(len(p["ids"]),size=n,replace=False,p=w)
                chosen.extend(p["ids"][idx].tolist()); chosen_x.append(p["X"][idx,:])
            key=tuple(sorted(map(str,chosen), key=lambda x:(not x.isdigit(),x)))
            if key in seen: duprej+=1; continue
            cx=np.vstack(chosen_x); av=abs_smd_from_moments(tmean,tvar,cx)
            if not bool(np.all(np.isfinite(av))) or bool(np.any(av>caliper)):
                balrej+=1; continue
            seen.add(key); vals.append(av); rep=len(seen)
            sf.write(f"{fam}\t{target}\t{rep}\t{','.join(key)}\n")
            mf.write(f"{fam}\t{target}\t{rep}\t{av[0]:.12g}\t{av[1]:.12g}\t{av[2]:.12g}\n")
            for g in key: reuse[(fam,target,g)]+=1
        if len(seen)!=B:
            raise SystemExit(f"HOLD: {fam}/{target} generated {len(seen)}/{B} sets after {attempts}; balance_rej={balrej}; dup={duprej}")
        arr=np.asarray(vals,float)
        for j,nm in enumerate(["log_gene_length","GC","annotation_degree"]):
            med=float(np.median(arr[:,j])); p95=float(np.quantile(arr[:,j],0.95))
            qc_rows.append({"family":fam,"target_id":target,"covariate":nm,
                            "median_abs_smd":med,"p95_abs_smd":p95,
                            "median_gate":med_gate,"p95_gate":p95_gate,
                            "status":"PASS" if med<=med_gate and p95<=p95_gate else "FAIL"})
        summary_rows.append({"family":fam,"target_id":target,"frozen_n":frozen_n,
                             "target_common_scored_n":len(target_ids),"common_scored_coverage":coverage,
                             "matched_set_n":B,"unique_set_n":len(seen),"generation_attempts":attempts,
                             "balance_rejections":balrej,"duplicate_rejections":duprej,
                             "acceptance_rate":B/attempts,"control_pool_n":len(control),
                             "global_hypothesis_union_excluded_n":len(all_hyp),
                             "matching_method":mc["matching_method"],"exact_chromosome_composition":"YES"})

qdf=pd.DataFrame(qc_rows)
if (qdf.status!="PASS").any(): raise SystemExit("HOLD: matched SMD quality gate failed")
pd.DataFrame(pool_rows).to_csv(OUT/"matched"/"MAR03_candidate_pool_qc.tsv",sep="\t",index=False)
pd.DataFrame(tilt_rows).to_csv(OUT/"matched"/"MAR03_matching_tilt_diagnostics.tsv",sep="\t",index=False)
qdf.to_csv(OUT/"MAR03_matching_qc.tsv",sep="\t",index=False)
pd.DataFrame(summary_rows).to_csv(OUT/"MAR03_matching_summary.tsv",sep="\t",index=False)
pd.DataFrame([{"family":f,"target_id":t,"control_gene_id":g,"selected_count":c,"selected_prop":c/B}
              for (f,t,g),c in sorted(reuse.items())]).to_csv(
              OUT/"matched"/"MAR03_matched_control_reuse.tsv",sep="\t",index=False)

# Matching is frozen before scores are loaded.
(OUT/"state"/"MAR03_MATCHING_PASS.txt").write_text(
    f"status=PASS\nstage=MAR03_MATCHING\nmethod={mc['matching_method']}\n"
    f"targets={len(targets)}\nB_per_target={B}\nscore_values_read=NO\nsame_matched_sets_PGC_SPARK=YES\n",
    encoding="utf-8")

# Now read gene scores and score the exact same matched sets in both datasets.
pgc=read_gene_scores(resolved["PGC_gene_results"]); spark=read_gene_scores(resolved["SPARK_gene_results"])
scores={"PGC":pgc,"SPARK":spark}

# Re-read matched sets into memory by target; 200k rows total is manageable.
matched=defaultdict(list)
with gzip.open(matched_path,"rt",encoding="utf-8") as f:
    r=csv.DictReader(f,delimiter="\t")
    for z in r:
        matched[(z["family"],z["target_id"])].append(tuple(z["genes"].split(",")))

res=[]
for fam,target in targets:
    target_ids=sorted(sets[(fam,target)] & set(d.index))
    sets_null=matched[(fam,target)]
    for ds,score in scores.items():
        usable=[g for g in target_ids if g in score]
        if len(usable)<2: raise SystemExit(f"HOLD: {ds} {fam}/{target} <2 scored target genes")
        obs=np.array([score[g] for g in usable],float)
        null_mean=[]; null_med=[]
        for s in sets_null:
            vals=[score[g] for g in s if g in score]
            if len(vals)!=len(s):
                raise SystemExit(f"HOLD: {ds} missing matched control scores for {fam}/{target}")
            null_mean.append(float(np.mean(vals))); null_med.append(float(np.median(vals)))
        nm=np.array(null_mean); nd=np.array(null_med)
        om=float(np.mean(obs)); od=float(np.median(obs))
        pmean=(1+int(np.sum(nm>=om)))/(B+1)
        pmed=(1+int(np.sum(nd>=od)))/(B+1)
        res.append({"dataset":ds,"family":fam,"target_id":target,
                    "frozen_n":len(sets[(fam,target)]),"common_scored_n":len(usable),
                    "observed_meanZ":om,"meanZ_percentile":float(np.mean(nm<=om)),"P_emp_meanZ":pmean,"q_BH_meanZ":"",
                    "observed_medianZ":od,"medianZ_percentile":float(np.mean(nd<=od)),"P_emp_medianZ":pmed,"q_BH_medianZ":"",
                    "matched_set_n":B,"meanZ_status":"PRIMARY_COMPLEMENTARY","medianZ_status":"SUPPORTIVE"})

# BH separately within dataset × family × statistic.
for ds in ["PGC","SPARK"]:
    for fam in ["A","B"]:
        idx=[i for i,z in enumerate(res) if z["dataset"]==ds and z["family"]==fam]
        q1=bh([res[i]["P_emp_meanZ"] for i in idx]); q2=bh([res[i]["P_emp_medianZ"] for i in idx])
        for i,q in zip(idx,q1): res[i]["q_BH_meanZ"]=q
        for i,q in zip(idx,q2): res[i]["q_BH_medianZ"]=q

write_tsv(OUT/"MAR03_matched_meanZ.tsv",
          ["dataset","family","target_id","frozen_n","common_scored_n",
           "observed_meanZ","meanZ_percentile","P_emp_meanZ","q_BH_meanZ",
           "observed_medianZ","medianZ_percentile","P_emp_medianZ","q_BH_medianZ",
           "matched_set_n","meanZ_status","medianZ_status"],res)
print("MAR03_MATCHED=PASS")
