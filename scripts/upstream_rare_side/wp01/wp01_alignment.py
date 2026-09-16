#!/usr/bin/env python3
from pathlib import Path
import hashlib, json, math, sys, numpy as np, pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from wp01_common import *
from wp01_lineage import verify_stable

VERSION="WP01_ALIGNMENT_STABLE_v1.2.0"
EXPECTED_GLOBAL_SCRIPT_SHA="42ba8c2dc34341925ec3c1a09b57dec1f72cf0ce5bfcf285c28aebf46a804b0c"
EXPECTED_GLOBAL_FILES={
  "WP01_global_overlap_summary.tsv":"39c979529f99b151e283aff60fdb1b1638f9eb2496135d4782ab41e4b3f7cb59",
  "WP01_matched_catalogue_balance.tsv":"3143947e5b19536aa6ed6f6c7bce29b7cdd956245d56cdeb91c109ff58f9e59d",
  "WP01_matched_catalogue_null.tsv.gz":"efeb38d33269cd75e2cc708230f2c3f215fac6ae8c56f3c4022eca3a89680e32",
  "WP01_global_match_tilt_fit.tsv":"b168311c112d585b9861f7452bde69605b8b5485ee066a86ad2d5148f4e01ba1",
  "WP01_global_match_weight_support.tsv.gz":"cab124d2121b693978499785eed182186952798376bbff823066fc22b854bc4f",
  "WP01_global_match_engine_diagnostics.tsv":"2f3c0263d4cd007249a5bbd1a15faa1926cbfb3eb54107921155c34c0f76aeda"
}

def parse_kv2(path):
    d={}
    for line in Path(path).read_text(errors="replace").splitlines():
        if "=" in line and not line.lstrip().startswith("#"):
            k,v=line.split("=",1); d[k.strip()]=v.strip()
    return d

def require_alignment_release(ma,out,locksha):
    # v1.2.0 uses one active implementation manifest. Historical v1.1.0/v1.1.1
    # release tokens remain immutable lineage evidence but are not compared
    # directly with the current script hash.
    return verify_stable(ma,out,locksha)

def verify_global_checkpoint(root,ma,out):
    rows=[]
    def add(k,status,details=""):
        rows.append({"check_id":k,"status":status,"details":details})

    # The global v1.0.2 code and lineage must be exactly the reviewed state.
    gscript=ma/"workflow/scripts/wp01_global_overlap.py"
    add("GLOBAL_CODE::v1.0.2",
        "PASS" if sha256_file(gscript)==EXPECTED_GLOBAL_SCRIPT_SHA else "FAIL",
        sha256_file(gscript))
    for rel in ["WP01_GLOBAL_MATCH_AMENDMENT_RELEASE.txt","WP01_GLOBAL_MATCH_HOTFIX_RELEASE.txt"]:
        p=out/rel
        add("GLOBAL_LINEAGE::"+rel,"PASS" if p.exists() and "decision=PASS" in p.read_text(errors="replace") else "FAIL",str(p))

    for f,h in EXPECTED_GLOBAL_FILES.items():
        p=out/f
        obs=sha256_file(p) if p.exists() else "MISSING"
        add("GLOBAL_HASH::"+f,"PASS" if obs==h else "FAIL",f"observed={obs}; expected={h}")

    try:
        s=pd.read_csv(out/"WP01_global_overlap_summary.tsv",sep="\t")
        row=s.iloc[0]
        count_ok=(int(row.N_universe)==19267 and int(row.IUIS_n)==500 and int(row.SFARI_R0_n)==935
                  and int(row.observed_overlap_n)==26 and int(row.matched_B)==10000)
        add("GLOBAL_SEMANTIC::counts","PASS" if count_ok else "FAIL",s.to_json())
        add("GLOBAL_SEMANTIC::balance_gate","PASS" if str(row.matched_balance_gate)=="PASS" else "FAIL",str(row.matched_balance_gate))
        # Recompute the two reported P values.
        ph=hypergeom_sf(19267,500,935,26)
        add("GLOBAL_SEMANTIC::hypergeom_recompute","PASS" if abs(ph-float(row.hypergeom_one_sided_P))<5e-12 else "FAIL",
            f"recomputed={ph}; reported={row.hypergeom_one_sided_P}")
    except Exception as e:
        add("GLOBAL_SEMANTIC::counts","FAIL",repr(e))

    try:
        null=pd.read_csv(out/"WP01_matched_catalogue_null.tsv.gz",sep="\t",compression="gzip")
        emp=(1+(null.overlap_n>=26).sum())/(len(null)+1)
        s=pd.read_csv(out/"WP01_global_overlap_summary.tsv",sep="\t").iloc[0]
        add("GLOBAL_NULL::rows","PASS" if len(null)==10000 else "FAIL",f"n={len(null)}")
        add("GLOBAL_NULL::unique_IUIS","PASS" if null.IUIS_set_sha256.nunique()==10000 else "FAIL",
            f"unique={null.IUIS_set_sha256.nunique()}")
        add("GLOBAL_NULL::unique_SFARI","PASS" if null.SFARI_set_sha256.nunique()==10000 else "FAIL",
            f"unique={null.SFARI_set_sha256.nunique()}")
        add("GLOBAL_NULL::empirical_P_recompute","PASS" if abs(emp-float(s.matched_empirical_P))<5e-12 else "FAIL",
            f"recomputed={emp}; reported={s.matched_empirical_P}")
        smdcols=[c for c in null.columns if c.endswith("_SMD")]
        maxs=max(float(null[c].abs().max()) for c in smdcols)
        add("GLOBAL_NULL::per_set_caliper","PASS" if maxs<=0.20+1e-12 else "FAIL",f"max_abs_SMD={maxs}")
    except Exception as e:
        add("GLOBAL_NULL::rows","FAIL",repr(e))

    try:
        bal=pd.read_csv(out/"WP01_matched_catalogue_balance.tsv",sep="\t")
        ok=(len(bal)==6 and (bal.status=="PASS").all()
            and (bal.median_abs_SMD<=0.10+1e-12).all()
            and (bal.p95_abs_SMD<=0.20+1e-12).all())
        add("GLOBAL_BALANCE::six_cells","PASS" if ok else "FAIL",bal.to_json())
    except Exception as e: add("GLOBAL_BALANCE::six_cells","FAIL",repr(e))

    try:
        d=pd.read_csv(out/"WP01_global_match_engine_diagnostics.tsv",sep="\t")
        ok=(len(d)==2 and (d.accepted_sets==10000).all() and (d.unique_set_hashes==10000).all()
            and np.isfinite(d.min_stratum_weight_ESS).all() and (d.min_stratum_weight_ESS>1).all()
            and (d.max_single_gene_weight<1).all())
        add("GLOBAL_ENGINE::diagnostics","PASS" if ok else "FAIL",d.to_json())
    except Exception as e: add("GLOBAL_ENGINE::diagnostics","FAIL",repr(e))

    try:
        t=pd.read_csv(out/"WP01_global_match_tilt_fit.tsv",sep="\t")
        ok=(len(t)==2 and t.newton_converged_1e8.map(truthy).all()
            and (t.max_abs_expected_mean_error_z<1e-8).all())
        add("GLOBAL_ENGINE::tilt_fit","PASS" if ok else "FAIL",t.to_json())
    except Exception as e: add("GLOBAL_ENGINE::tilt_fit","FAIL",repr(e))

    df=pd.DataFrame(rows)
    df.to_csv(out/"WP01_global_checkpoint_checks.tsv",sep="\t",index=False)
    bad=df[df.status!="PASS"]
    obj={
      "stage":"WP01","checkpoint":"GLOBAL_MATCH_v1.0.2","status":"PASS" if len(bad)==0 else "HOLD",
      "preanalysis_lock_sha256":sha256_file(out/"WP01_preanalysis_lock.json"),
      "global_script_sha256":sha256_file(gscript),
      "files":{f:sha256_file(out/f) for f in EXPECTED_GLOBAL_FILES},
      "checks_n":len(df),"failures_n":len(bad)
    }
    write_json(out/"WP01_GLOBAL_CHECKPOINT_v1.0.2.json",obj)
    marker=out/("WP01_GLOBAL_CHECKPOINT_PASS.txt" if len(bad)==0 else "WP01_GLOBAL_CHECKPOINT_HOLD.txt")
    marker.write_text(obj["status"]+"\n",encoding="utf-8")
    if len(bad):
        raise RuntimeError("Reviewed global v1.0.2 checkpoint failed revalidation")
    return obj

def build_canonical_frame(root,ma,out):
    rows=[]
    def add(k,obs,exp,details=""):
        status="PASS" if obs==exp else "FAIL"
        rows.append({"check_id":k,"observed":obs,"expected":exp,"status":status,"details":details})

    # Universe is the canonical row-level source. It already contains all flags/covariates.
    u=pd.read_csv(root/"01_gene_universe/01_gene_universe.tsv",sep="\t",low_memory=False)
    required=["HGNC_id","HGNC_symbol","chr","start","end","log_gene_length","GC_fraction","GC_z",
              "IUIS_flag","SFARI_R0_highconf_flag","SFARI_R1_highconf_flag","CoreSeed_flag","R0_group",
              "annotation_degree","annotation_degree_log1p","annotation_degree_log1p_z"]
    miss=[c for c in required if c not in u.columns]
    if miss: raise RuntimeError("01_gene_universe canonical columns missing: "+",".join(miss))
    add("universe_rows",len(u),19267)
    add("universe_unique_symbol",u.HGNC_symbol.nunique(),19267)
    add("universe_unique_HGNC_id",u.HGNC_id.nunique(),19267)

    # Independent gene-group cross-check. No merge-back into the canonical frame.
    g=pd.read_csv(root/"01_gene_universe/02_gene_groups.tsv",sep="\t",low_memory=False)
    greq=["HGNC_id","approved_symbol","IUIS_flag","SFARI_R0_highconf_flag","SFARI_R1_highconf_flag","CoreSeed_flag","R0_group"]
    miss=[c for c in greq if c not in g.columns]
    if miss: raise RuntimeError("02_gene_groups canonical columns missing: "+",".join(miss))
    gg=u[["HGNC_id","HGNC_symbol","IUIS_flag","SFARI_R0_highconf_flag","SFARI_R1_highconf_flag","CoreSeed_flag","R0_group"]].merge(
        g,left_on=["HGNC_id","HGNC_symbol"],right_on=["HGNC_id","approved_symbol"],
        suffixes=("_universe","_groups"),validate="one_to_one")
    add("gene_groups_rows",len(gg),19267)
    for c in ["IUIS_flag","SFARI_R0_highconf_flag","SFARI_R1_highconf_flag","CoreSeed_flag","R0_group"]:
        n=int((gg[f"{c}_universe"].astype(str)!=gg[f"{c}_groups"].astype(str)).sum())
        add("gene_group_mismatch::"+c,n,0)

    # Independent annotation table cross-check.
    a=pd.read_csv(root/"02_pathways/annotation_degree.tsv",sep="\t",low_memory=False)
    areq=["HGNC_id","approved_symbol","annotation_degree","annotation_degree_log1p","annotation_degree_log1p_z"]
    miss=[c for c in areq if c not in a.columns]
    if miss: raise RuntimeError("annotation_degree canonical columns missing: "+",".join(miss))
    aa=u[["HGNC_id","HGNC_symbol","annotation_degree","annotation_degree_log1p","annotation_degree_log1p_z"]].merge(
        a,left_on=["HGNC_id","HGNC_symbol"],right_on=["HGNC_id","approved_symbol"],
        suffixes=("_universe","_annotation"),validate="one_to_one")
    for c in ["annotation_degree","annotation_degree_log1p","annotation_degree_log1p_z"]:
        md=float(np.nanmax(np.abs(pd.to_numeric(aa[f"{c}_universe"])-pd.to_numeric(aa[f"{c}_annotation"]))))
        add("annotation_max_abs_diff::"+c,round(md,15),0.0)

    # Legacy model basis becomes the single authoritative MDV3 design interface.
    # IMPORTANT: normalize the join key explicitly before any alignment.
    # Do NOT rely on implicit merge suffixes or assume both tables already use "gene".
    sp=pd.read_csv(root/"04_bias_empirical/MDV3_spline_basis.tsv.gz",sep="\t",compression="gzip",low_memory=False)
    sreq=["HGNC_symbol","length_ns1","length_ns2","length_ns3","GC_z_model","annotation_z_model"]
    miss=[c for c in sreq if c not in sp.columns]
    if miss: raise RuntimeError("MDV3 spline basis columns missing: "+",".join(miss))
    sp=sp[sreq].copy().rename(columns={"HGNC_symbol":"gene"})
    sp["gene"]=sp["gene"].astype(str).str.strip()
    add("spline_rows",len(sp),19267)
    add("spline_unique_gene",sp.gene.nunique(),19267)
    add("spline_duplicate_gene_n",int(sp.gene.duplicated().sum()),0)

    f=u[required].copy().rename(columns={
        "HGNC_symbol":"gene","GC_fraction":"GC","annotation_degree_log1p":"log_annotation",
        "GC_z":"z_GC","annotation_degree_log1p_z":"z_log_annotation"
    })
    f["gene"]=f["gene"].astype(str).str.strip()

    # Key-set identity is checked before value alignment.
    f_gene_set=set(f["gene"])
    sp_gene_set=set(sp["gene"])
    add("spline_missing_from_universe_n",len(sp_gene_set-f_gene_set),0)
    add("universe_missing_from_spline_n",len(f_gene_set-sp_gene_set),0)
    if sp["gene"].duplicated().any():
        raise RuntimeError("MDV3 spline basis has duplicate gene keys")
    if f_gene_set!=sp_gene_set:
        raise RuntimeError(
            f"MDV3 spline/universe gene-key identity failure: "
            f"spline_only={len(sp_gene_set-f_gene_set)}, universe_only={len(f_gene_set-sp_gene_set)}"
        )

    # Align by explicit key reindexing rather than merge(). This preserves the
    # canonical universe row order and prevents .x/.y or missing-key surprises.
    basis_cols=["length_ns1","length_ns2","length_ns3","GC_z_model","annotation_z_model"]
    sp_idx=sp.set_index("gene",verify_integrity=True)
    aligned=sp_idx.loc[f["gene"],basis_cols]
    if len(aligned)!=len(f):
        raise RuntimeError("Legacy model basis alignment changed row count")
    if aligned.isna().any().any():
        raise RuntimeError("Missing legacy model basis after explicit key alignment")
    for c in basis_cols:
        f[c]=aligned[c].to_numpy()

    # Post-alignment identity guard: row order and model-vector lengths must match.
    add("spline_aligned_rows",len(aligned),19267)
    dgc=float(np.max(np.abs(f.z_GC-f.GC_z_model)))
    dan=float(np.max(np.abs(f.z_log_annotation-f.annotation_z_model)))
    rows.append({"check_id":"legacy_model_GC_z_max_abs_diff","observed":dgc,"expected":"<=1e-12",
                 "status":"PASS" if dgc<=1e-12 else "FAIL","details":""})
    rows.append({"check_id":"legacy_model_annotation_z_max_abs_diff","observed":dan,"expected":"<=1e-12",
                 "status":"PASS" if dan<=1e-12 else "FAIL","details":""})

    f["raw26_flag"]=((f.IUIS_flag.astype(int)==1)&(f.SFARI_R0_highconf_flag.astype(int)==1)).astype(int)
    f["z_log_gene_length"]=(f.log_gene_length-f.log_gene_length.mean())/f.log_gene_length.std(ddof=1)
    # R type=7 quantiles == numpy linear quantile.
    bounds=np.quantile(f.log_gene_length.to_numpy(dtype=float),[.2,.4,.6,.8],method="linear")
    f["length_quintile"]=np.searchsorted(bounds,f.log_gene_length.to_numpy(dtype=float),side="right")+1

    add("IUIS_n",int(f.IUIS_flag.sum()),500)
    add("SFARI_R0_n",int(f.SFARI_R0_highconf_flag.sum()),935)
    add("raw26_n",int(f.raw26_flag.sum()),26)
    add("CoreSeed25_n",int(f.CoreSeed_flag.sum()),25)
    add("CoreSeed_outside_raw26_n",int(((f.CoreSeed_flag==1)&(f.raw26_flag==0)).sum()),0)

    # Frozen raw26/Core25 identity cross-check.
    rawf=pd.read_csv(out/"WP01_raw26_gene_set.tsv",sep="\t")
    coref=pd.read_csv(out/"WP01_core25_historical_gene_set.tsv",sep="\t")
    add("raw26_set_identity",int(set(f.loc[f.raw26_flag==1,"gene"])==set(rawf.gene)),1)
    add("core25_set_identity",int(set(f.loc[f.CoreSeed_flag==1,"gene"])==set(coref.gene)),1)
    add("raw26_minus_core25_is_C4B",int(set(f.loc[f.raw26_flag==1,"gene"])-set(f.loc[f.CoreSeed_flag==1,"gene"])=={"C4B"}),1)

    # Pathway identity must agree across manifest, membership and legacy MDV3 results.
    pm=pd.read_csv(root/"02_pathways/pathway_manifest.tsv",sep="\t",low_memory=False)
    mem=pd.read_csv(root/"02_pathways/03_pathway_membership.tsv.gz",sep="\t",compression="gzip",low_memory=False)
    legacy=pd.read_csv(root/"04_bias_empirical/04_pathway_all.tsv.gz",sep="\t",compression="gzip",low_memory=False)
    prim=set(map(tuple,pm.loc[pm.primary_10_500_flag==1,["source","pathway_id"]].drop_duplicates().to_numpy()))
    mids=set(map(tuple,mem[["source","pathway_id"]].drop_duplicates().to_numpy()))
    lids=set(map(tuple,legacy[["source","pathway_id"]].drop_duplicates().to_numpy()))
    add("primary_pathway_n",len(prim),6671)
    add("GO_BP_n",sum(1 for s,p in prim if s=="GO_BP"),5450)
    add("Reactome_n",sum(1 for s,p in prim if s=="Reactome"),1221)
    add("manifest_membership_symdiff",len(prim.symmetric_difference(mids)),0)
    add("manifest_legacy_mdv3_symdiff",len(prim.symmetric_difference(lids)),0)
    add("membership_genes_missing_universe",len(set(mem.approved_symbol)-set(f.gene)),0)

    audit=pd.DataFrame(rows)
    audit.to_csv(out/"WP01_canonical_schema_audit.tsv",sep="\t",index=False)
    bad=audit[audit.status!="PASS"]
    if len(bad):
        raise RuntimeError("Canonical schema audit failed:\n"+bad.to_string(index=False))

    cols=["HGNC_id","gene","chr","start","end","IUIS_flag","SFARI_R0_highconf_flag",
          "SFARI_R1_highconf_flag","CoreSeed_flag","raw26_flag","R0_group",
          "log_gene_length","GC","annotation_degree","log_annotation",
          "z_log_gene_length","z_GC","z_log_annotation","length_quintile",
          "length_ns1","length_ns2","length_ns3","GC_z_model","annotation_z_model"]
    f[cols].to_csv(out/"WP01_model_frame.tsv.gz",sep="\t",index=False,compression="gzip")
    frame_sha=sha256_file(out/"WP01_model_frame.tsv.gz")
    (out/"WP01_model_frame.sha256").write_text(frame_sha+"  WP01_model_frame.tsv.gz\n",encoding="utf-8")

    schema={
      "version":VERSION,"rows":len(f),"columns":cols,
      "primary_anchor":"raw26_flag","historical_sensitivity":"CoreSeed_flag",
      "matching_features":["z_log_gene_length","z_GC","z_log_annotation","length_quintile"],
      "mdv3_model_features":["length_ns1","length_ns2","length_ns3","GC_z_model","annotation_z_model"],
      "model_frame_sha256":frame_sha,
      "length_quintile_bounds":bounds.tolist()
    }
    write_json(out/"WP01_canonical_schema.json",schema)
    return frame_sha

def main():
    if len(sys.argv)<2: raise SystemExit("usage: wp01_alignment.py CH3_ROOT")
    root=Path(sys.argv[1]).resolve()
    ma=root/"12_molecular_autism_revision"
    out=ma/"01_anchor_ingest/WP01_raw26_primary_v1"
    locksha=sha256_file(out/"WP01_preanalysis_lock.json")
    require_alignment_release(ma,out,locksha)

    global_obj=verify_global_checkpoint(root,ma,out)
    frame_sha=build_canonical_frame(root,ma,out)

    # Freeze the aligned code and checkpoints after successful construction.
    code_files=[
      ma/"workflow/scripts/wp01_common.py",
      ma/"workflow/scripts/wp01_lineage.py",
      ma/"workflow/scripts/wp01_alignment.py",
      ma/"workflow/scripts/wp01_raw26_mdv3.R",
      ma/"workflow/scripts/wp01_matching.py",
      ma/"workflow/scripts/wp01_stage1_validate.py",
      ma/"workflow/scripts/wp01_empirical_tiering.py",
      ma/"workflow/scripts/wp01_delta_core25.py",
      ma/"workflow/scripts/wp01_postflight.py",
      ma/"bin/run_wp01.sh",ma/"bin/build_wp01_review_bundle.sh",ma/"bin/wp01_release_gate.sh",
      ma/"config/wp01_contract.json"
    ]
    code={str(p.relative_to(ma)):sha256_file(p) for p in code_files}
    write_json(out/"WP01_ALIGNMENT_CODE_LOCK.json",{"version":VERSION,"files":code})
    alignment_lock={
      "stage":"WP01","version":VERSION,"status":"PASS_CANDIDATE",
      "preanalysis_lock_sha256":locksha,
      "global_checkpoint_sha256":sha256_file(out/"WP01_GLOBAL_CHECKPOINT_v1.0.2.json"),
      "model_frame_sha256":frame_sha,
      "canonical_schema_sha256":sha256_file(out/"WP01_canonical_schema.json"),
      "canonical_schema_audit_sha256":sha256_file(out/"WP01_canonical_schema_audit.tsv"),
      "alignment_code_lock_sha256":sha256_file(out/"WP01_ALIGNMENT_CODE_LOCK.json")
    }
    write_json(out/"WP01_ALIGNMENT_LOCK.json",alignment_lock)
    (out/"WP01_ALIGNMENT_PASS.txt").write_text("PASS\n",encoding="utf-8")
    print(json.dumps(alignment_lock,indent=2))

if __name__=="__main__":
    main()
