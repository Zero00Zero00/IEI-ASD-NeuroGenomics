#!/usr/bin/env python3
from pathlib import Path
import json,sys,math,hashlib,numpy as np,pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from wp01_common import *

ENGINE_VERSION = "global-match-v1.0.2-stratified-exponential-tilt-interface-hotfix"

def qbin(s,q):
    # Preserve the v1.0.0 quantile-stratum definition.
    return pd.qcut(s.rank(method="first"),q=q,labels=False,duplicates="drop").astype(int)

def parse_kv_local(path: Path):
    d={}
    for line in path.read_text(errors="replace").splitlines():
        if "=" in line and not line.lstrip().startswith("#"):
            k,v=line.split("=",1); d[k.strip()]=v.strip()
    return d

def require_global_match_amendment(out: Path, locksha: str):
    # Preserve the approved v1.0.1 matching-engine amendment as the parent decision.
    rel1=out/"WP01_GLOBAL_MATCH_AMENDMENT_RELEASE.txt"
    if not rel1.exists():
        raise RuntimeError("v1.0.1 global-match amendment release missing")
    r1=parse_kv_local(rel1)
    if r1.get("decision")!="PASS" or r1.get("stage")!="WP01" or r1.get("amendment")!="GLOBAL_MATCH_ENGINE_v1.0.1":
        raise RuntimeError("Invalid v1.0.1 global-match amendment release token")
    if r1.get("preanalysis_lock_sha256")!=locksha:
        raise RuntimeError("v1.0.1 amendment bound to a different preanalysis lock")
    if r1.get("patched_global_script_sha256")!="1439eda02941de3f5b07ca9236e412c86917e65b5925300f64086dda94aec101":
        raise RuntimeError("v1.0.1 amendment token does not bind the approved parent script")

    # v1.0.2 is a pure Python interface hotfix and requires a second human release.
    ready2=out/"WP01_GLOBAL_MATCH_HOTFIX_READY.txt"
    rel2=out/"WP01_GLOBAL_MATCH_HOTFIX_RELEASE.txt"
    if not ready2.exists() or not rel2.exists():
        raise RuntimeError("v1.0.2 global-match hotfix READY/release missing")
    current_script_sha=sha256_file(Path(__file__))
    r2=parse_kv_local(rel2)
    if r2.get("decision")!="PASS" or r2.get("stage")!="WP01" or r2.get("amendment")!="GLOBAL_MATCH_INTERFACE_HOTFIX_v1.0.2":
        raise RuntimeError("Invalid v1.0.2 global-match hotfix release token")
    if r2.get("preanalysis_lock_sha256")!=locksha:
        raise RuntimeError("v1.0.2 hotfix bound to a different preanalysis lock")
    if r2.get("parent_global_script_sha256")!="1439eda02941de3f5b07ca9236e412c86917e65b5925300f64086dda94aec101":
        raise RuntimeError("v1.0.2 hotfix parent-script hash mismatch")
    if r2.get("patched_global_script_sha256")!=current_script_sha:
        raise RuntimeError("v1.0.2 hotfix bound to a different current script")

def _strata_objects(frame, mask, strata_col, feat_cols):
    obs=frame.loc[mask]
    counts=obs.groupby(strata_col).size().to_dict()
    objects=[]
    for st,n in sorted(counts.items(),key=lambda x:str(x[0])):
        cand=frame.loc[frame[strata_col]==st]
        Z=cand[feat_cols].to_numpy(dtype=float)
        idx=cand.index.to_numpy(dtype=int)
        if n>len(idx):
            raise RuntimeError(f"Stratum {st}: need {n}, have {len(idx)}")
        objects.append((st,int(n),idx,Z))
    return obs,counts,objects

def _expected_moments(lam, objects, n_total):
    mu=np.zeros(3,dtype=float)
    cov=np.zeros((3,3),dtype=float)
    weight_diag=[]
    for st,n,idx,Z in objects:
        eta=Z@lam
        eta=eta-np.max(eta)
        w=np.exp(np.clip(eta,-60,60))
        p=w/w.sum()
        mus=p@Z
        C=Z-mus
        covs=(C*p[:,None]).T@C
        frac=n/n_total
        mu += frac*mus
        cov += frac*covs
        ess=1.0/np.sum(p*p)
        weight_diag.append((st,n,len(idx),ess,float(np.max(p))))
    return mu,cov,weight_diag

def fit_stratified_tilt(frame, mask, strata_col, feat_cols, label):
    obs,counts,objects=_strata_objects(frame,mask,strata_col,feat_cols)
    target=obs[feat_cols].mean().to_numpy(dtype=float)
    lam=np.zeros(len(feat_cols),dtype=float)
    converged=False
    history=[]
    for it in range(1,101):
        mu,cov,wdiag=_expected_moments(lam,objects,len(obs))
        diff=target-mu
        err=float(np.max(np.abs(diff)))
        history.append((it,err,*lam.tolist()))
        if err<1e-8:
            converged=True; break
        ridge=1e-6*np.eye(len(feat_cols))
        try:
            delta=np.linalg.solve(cov+ridge,diff)
        except np.linalg.LinAlgError:
            delta=np.linalg.pinv(cov+ridge)@diff
        dn=np.linalg.norm(delta)
        if dn>5:
            delta=delta*(5/dn)
        olderr=err
        accepted=False
        for step in [1.0,0.5,0.25,0.1,0.05,0.01]:
            trial=lam+step*delta
            tmu,_,_=_expected_moments(trial,objects,len(obs))
            terr=float(np.max(np.abs(target-tmu)))
            if np.isfinite(terr) and terr<olderr:
                lam=trial; accepted=True; break
        if not accepted:
            break
    mu,cov,wdiag=_expected_moments(lam,objects,len(obs))
    final_err=float(np.max(np.abs(target-mu)))
    if not np.all(np.isfinite(lam)) or not np.isfinite(final_err):
        raise RuntimeError(f"{label}: non-finite exponential-tilt solution")
    # Sampling can still be valid when Newton stops slightly above numerical tolerance;
    # structural acceptance below, not the solver tolerance, is the inferential gate.
    diag=pd.DataFrame(wdiag,columns=["stratum","target_n","candidate_n","weight_ESS","max_single_gene_weight"])
    diag["catalogue"]=label
    fit=pd.DataFrame([{
        "catalogue":label,
        "lambda_log_gene_length":lam[0],
        "lambda_GC":lam[1],
        "lambda_log_annotation":lam[2],
        "target_mean_z_log_gene_length":target[0],
        "target_mean_GC_z":target[1],
        "target_mean_anno_z":target[2],
        "fitted_expected_mean_z_log_gene_length":mu[0],
        "fitted_expected_mean_GC_z":mu[1],
        "fitted_expected_mean_anno_z":mu[2],
        "max_abs_expected_mean_error_z":final_err,
        "newton_converged_1e8":converged,
        "iterations":len(history)
    }])
    return obs,objects,lam,diag,fit

def make_probabilities(lam,objects):
    out={}
    for st,n,idx,Z in objects:
        eta=Z@lam; eta=eta-np.max(eta)
        w=np.exp(np.clip(eta,-60,60)); p=w/w.sum()
        out[st]=(n,idx,p)
    return out

def draw_one_catalogue(rng, frame, probs, obs, raw_cols, per_set_caliper, seen):
    tries=0
    while True:
        tries+=1
        picked=[]
        for st,(n,idx,p) in probs.items():
            picked.extend(rng.choice(idx,size=n,replace=False,p=p).tolist())
        picked=np.array(picked,dtype=int)
        h=hashlib.sha256(("\n".join(map(str,sorted(picked.tolist())))+"\n").encode()).hexdigest()
        if h in seen:
            continue
        smds={}
        for col in raw_cols:
            smds[col]=smd_pooled(frame.loc[picked,col],obs[col])
        maxabs=max(abs(v) for v in smds.values())
        if np.isfinite(maxabs) and maxabs<=per_set_caliper:
            seen.add(h)
            return picked,h,smds,tries
        if tries>5000:
            raise RuntimeError("Unable to generate a structurally accepted catalogue within 5000 attempts")

def main():
    root=Path(sys.argv[1]).resolve(); ma=root/"12_molecular_autism_revision"
    out=ma/"01_anchor_ingest/WP01_raw26_primary_v1"
    locksha=sha256_file(out/"WP01_preanalysis_lock.json")
    require_release(out/"WP01_STOPA_RELEASE.txt",locksha,"A")
    require_global_match_amendment(out,locksha)

    c=read_contract(ma)
    B=int(c["rules"]["global_matched_null_B"])
    q=int(c["rules"]["global_joint_quantile_bins"])
    seed=int(c["seeds"]["global_catalogue_null"])
    median_gate=float(c["rules"]["global_balance_median_gate"])
    p95_gate=float(c["rules"]["global_balance_p95_gate"])
    # Added as a QC-only acceptance rule; it is equal to the already frozen P95 gate
    # and therefore does not relax any prespecified balance threshold.
    per_set_caliper=p95_gate

    frame,meta=load_structural_frame(root)
    frame=frame.copy()
    for col in ["log_gene_length","GC","log_annotation"]:
        frame[col+"_bin"]=qbin(frame[col],q)
    frame["joint_stratum"]=(
        frame["log_gene_length_bin"].astype(str)+"|"+
        frame["GC_bin"].astype(str)+"|"+
        frame["log_annotation_bin"].astype(str)
    )

    iuis=frame.IUIS_flag.to_numpy()
    sfari=frame.SFARI_R0_highconf_flag.to_numpy()
    N=len(frame); K=int(iuis.sum()); n=int(sfari.sum()); x=int((iuis & sfari).sum())
    expected=K*n/N; p_h=hypergeom_sf(N,K,n,x)

    # load_structural_frame() creates Python-standardized columns as:
    # z_log_gene_length, z_GC, z_log_annotation.
    # v1.0.1 incorrectly used the R-side aliases GC_z and anno_z.
    feat_cols=["z_log_gene_length","z_GC","z_log_annotation"]
    raw_cols=["log_gene_length","GC","log_annotation"]
    missing_feat=[c for c in feat_cols if c not in frame.columns]
    if missing_feat:
        raise RuntimeError(
            "Standardized structural feature interface unresolved: "
            + ",".join(missing_feat)
            + "; available standardized columns="
            + ",".join([c for c in frame.columns if c.startswith("z_")])
        )

    obs_i,obj_i,lam_i,wdiag_i,fit_i=fit_stratified_tilt(frame,iuis,"joint_stratum",feat_cols,"IUIS")
    obs_s,obj_s,lam_s,wdiag_s,fit_s=fit_stratified_tilt(frame,sfari,"joint_stratum",feat_cols,"SFARI")
    probs_i=make_probabilities(lam_i,obj_i)
    probs_s=make_probabilities(lam_s,obj_s)

    # Independent deterministic RNG streams from the unchanged frozen seed.
    ss=np.random.SeedSequence(seed)
    rng_i,rng_s=[np.random.default_rng(z) for z in ss.spawn(2)]
    seen_i=set(); seen_s=set()
    rec=[]
    total_try_i=0; total_try_s=0

    for b in range(1,B+1):
        I,hash_i,smd_i,tries_i=draw_one_catalogue(rng_i,frame,probs_i,obs_i,raw_cols,per_set_caliper,seen_i)
        S,hash_s,smd_s,tries_s=draw_one_catalogue(rng_s,frame,probs_s,obs_s,raw_cols,per_set_caliper,seen_s)
        total_try_i+=tries_i; total_try_s+=tries_s
        # Overlap is scored only after both catalogues passed the structural caliper.
        overlap=len(set(I.tolist()).intersection(S.tolist()))
        r={
            "replicate_id":b,
            "IUIS_set_sha256":hash_i,
            "SFARI_set_sha256":hash_s,
            "IUIS_draw_attempts":tries_i,
            "SFARI_draw_attempts":tries_s,
            "overlap_n":overlap
        }
        for col in raw_cols:
            r[f"IUIS_{col}_SMD"]=smd_i[col]
            r[f"SFARI_{col}_SMD"]=smd_s[col]
        rec.append(r)
        if b%1000==0:
            print(f"{ENGINE_VERSION}: accepted {b}/{B}",flush=True)

    null=pd.DataFrame(rec)
    if null.IUIS_set_sha256.nunique()!=B or null.SFARI_set_sha256.nunique()!=B:
        raise RuntimeError("Accepted pseudo-catalogue sets are not unique")

    balance=[]
    for label in ["IUIS","SFARI"]:
        for col in raw_cols:
            a=np.abs(null[f"{label}_{col}_SMD"].to_numpy(dtype=float))
            balance.append({
                "catalogue":label,"covariate":col,
                "median_abs_SMD":float(np.median(a)),
                "p95_abs_SMD":float(np.quantile(a,.95)),
                "max_abs_SMD":float(np.max(a)),
                "median_gate":median_gate,
                "p95_gate":p95_gate,
                "status":"PASS" if (np.median(a)<=median_gate and np.quantile(a,.95)<=p95_gate) else "FAIL"
            })
    bal=pd.DataFrame(balance)

    # Structural gate is evaluated before the overlap empirical P is released.
    gate=(bal.status=="PASS").all()
    null.to_csv(out/"WP01_matched_catalogue_null.tsv.gz",sep="\t",index=False,compression="gzip")
    bal.to_csv(out/"WP01_matched_catalogue_balance.tsv",sep="\t",index=False)

    fit=pd.concat([fit_i,fit_s],ignore_index=True)
    fit["engine_version"]=ENGINE_VERSION
    fit["joint_bins_per_covariate"]=q
    fit["per_set_abs_SMD_caliper"]=per_set_caliper
    fit.to_csv(out/"WP01_global_match_tilt_fit.tsv",sep="\t",index=False)

    wdiag=pd.concat([wdiag_i,wdiag_s],ignore_index=True)
    wdiag["engine_version"]=ENGINE_VERSION
    wdiag.to_csv(out/"WP01_global_match_weight_support.tsv.gz",sep="\t",index=False,compression="gzip")

    diagnostics=pd.DataFrame([
        {"catalogue":"IUIS","accepted_sets":B,"total_draw_attempts":total_try_i,
         "acceptance_rate":B/total_try_i,"unique_set_hashes":null.IUIS_set_sha256.nunique(),
         "min_stratum_weight_ESS":float(wdiag_i.weight_ESS.min()),
         "median_stratum_weight_ESS":float(wdiag_i.weight_ESS.median()),
         "max_single_gene_weight":float(wdiag_i.max_single_gene_weight.max())},
        {"catalogue":"SFARI","accepted_sets":B,"total_draw_attempts":total_try_s,
         "acceptance_rate":B/total_try_s,"unique_set_hashes":null.SFARI_set_sha256.nunique(),
         "min_stratum_weight_ESS":float(wdiag_s.weight_ESS.min()),
         "median_stratum_weight_ESS":float(wdiag_s.weight_ESS.median()),
         "max_single_gene_weight":float(wdiag_s.max_single_gene_weight.max())}
    ])
    diagnostics["engine_version"]=ENGINE_VERSION
    diagnostics["per_set_abs_SMD_caliper"]=per_set_caliper
    diagnostics.to_csv(out/"WP01_global_match_engine_diagnostics.tsv",sep="\t",index=False)

    if not gate:
        # Do not calculate/release matched empirical overlap P when the structural null fails.
        hold=pd.DataFrame([{
            "N_universe":N,"IUIS_n":K,"SFARI_R0_n":n,"observed_overlap_n":x,
            "hypergeom_expected_overlap":expected,"hypergeom_fold_enrichment":x/expected,
            "hypergeom_one_sided_P":p_h,"matched_B":B,
            "matched_empirical_P":np.nan,"matched_balance_gate":"FAIL",
            "seed":seed,"joint_bins_per_covariate":q,"engine_version":ENGINE_VERSION
        }])
        hold.to_csv(out/"WP01_global_overlap_summary.tsv",sep="\t",index=False)
        raise SystemExit("Matched-catalogue balance gate failed after v1.0.1 amendment")

    # Only after structural QC passes do we score the matched overlap null.
    p_emp=(1+(null.overlap_n>=x).sum())/(B+1)
    summary=pd.DataFrame([{
        "N_universe":N,"IUIS_n":K,"SFARI_R0_n":n,"observed_overlap_n":x,
        "hypergeom_expected_overlap":expected,"hypergeom_fold_enrichment":x/expected,
        "hypergeom_one_sided_P":p_h,"matched_B":B,
        "matched_overlap_mean":null.overlap_n.mean(),
        "matched_overlap_median":null.overlap_n.median(),
        "matched_overlap_p95":null.overlap_n.quantile(.95),
        "matched_empirical_P":p_emp,
        "matched_balance_gate":"PASS","seed":seed,
        "joint_bins_per_covariate":q,
        "engine_version":ENGINE_VERSION,
        "per_set_abs_SMD_caliper":per_set_caliper,
        "IUIS_acceptance_rate":B/total_try_i,
        "SFARI_acceptance_rate":B/total_try_s
    }])
    summary.to_csv(out/"WP01_global_overlap_summary.tsv",sep="\t",index=False)
    print(summary.to_string(index=False))

if __name__=="__main__":
    main()
