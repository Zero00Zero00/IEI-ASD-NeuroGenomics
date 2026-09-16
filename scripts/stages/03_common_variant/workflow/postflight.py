
from common import *
import shutil

ensure_dirs()
checks=[]
def ck(cid,ok,obs,exp="",detail=""):
    checks.append({"check_id":cid,"status":"PASS" if ok else "HOLD","observed":obs,"expected":exp,"detail":detail})

req=[
 "MAR03_gene_set_results_all.tsv","MAR03_pgc_higher_order.tsv","MAR03_spark_higher_order.tsv",
 "MAR03_standalone_tests.tsv","MAR03_matched_meanZ.tsv","MAR03_matching_qc.tsv",
 "MAR03_matching_summary.tsv","MAR03_sensitivity.tsv","MAR03_preanalysis_coverage.tsv",
 "MAR03_preanalysis_lock.json","MAR03_STOPA_RELEASE.txt"
]
for f in req: ck("exists_"+f,(OUT/f).is_file(),f,"exists")

if all((OUT/f).is_file() for f in req):
    _,allr=read_tsv(OUT/"MAR03_gene_set_results_all.tsv")
    ck("PGC_A_marginal_n",sum(z["dataset"]=="PGC" and z["family"]=="A" and z["model"]=="MARGINAL" for z in allr)==7,
       sum(z["dataset"]=="PGC" and z["family"]=="A" and z["model"]=="MARGINAL" for z in allr),7)
    ck("SPARK_A_marginal_n",sum(z["dataset"]=="SPARK" and z["family"]=="A" and z["model"]=="MARGINAL" for z in allr)==7,
       sum(z["dataset"]=="SPARK" and z["family"]=="A" and z["model"]=="MARGINAL" for z in allr),7)
    ck("standalone_marginal_n",sum(z["family"]=="B" and z["model"]=="MARGINAL" for z in allr)==26,
       sum(z["family"]=="B" and z["model"]=="MARGINAL" for z in allr),26)
    ck("conditional_A_n",sum(z["family"]=="A" and z["model"].startswith("CONDITIONAL") for z in allr)==14,
       sum(z["family"]=="A" and z["model"].startswith("CONDITIONAL") for z in allr),14)

    _,mr=read_tsv(OUT/"MAR03_matched_meanZ.tsv")
    ck("matched_result_n",len(mr)==40,len(mr),40)
    ck("matched_B",all(int(z["matched_set_n"])==10000 for z in mr),
       min(int(z["matched_set_n"]) for z in mr),10000)
    _,ms=read_tsv(OUT/"MAR03_matching_summary.tsv")
    ck("matched_target_n",len(ms)==20,len(ms),20)
    ck("matched_unique_sets",all(int(z["unique_set_n"])==10000 for z in ms),
       min(int(z["unique_set_n"]) for z in ms),10000)
    _,mq=read_tsv(OUT/"MAR03_matching_qc.tsv")
    ck("matching_quality_all_PASS",all(z["status"]=="PASS" for z in mq),
       sum(z["status"]=="PASS" for z in mq),len(mq))

    _,sr=read_tsv(OUT/"MAR03_sensitivity.tsv")
    ck("sensitivity_n",len(sr)==70,len(sr),70)

    _,cr=read_tsv(OUT/"MAR03_preanalysis_coverage.tsv")
    ck("dataset_mapping_nonzero",all(int(z["mapped_n"])>0 for z in cr),
       min(int(z["mapped_n"]) for z in cr),">0")

    tok={}
    for line in (OUT/"MAR03_STOPA_RELEASE.txt").read_text().splitlines():
        if "=" in line:
            k,v=line.split("=",1); tok[k]=v
    ck("release_lock_match",tok.get("preanalysis_lock_sha256")==sha256_file(OUT/"MAR03_preanalysis_lock.json"),
       tok.get("preanalysis_lock_sha256",""),sha256_file(OUT/"MAR03_preanalysis_lock.json"))

write_tsv(OUT/"MAR03_postflight_checks.tsv",["check_id","status","observed","expected","detail"],checks)
bad=[z for z in checks if z["status"]!="PASS"]
if bad:
    print(f"HOLD: {len(bad)} postflight checks failed. See {OUT/'MAR03_postflight_checks.tsv'}")
    raise SystemExit(40)

# Canonical source-data.
_,labels=read_tsv(PKG/"config"/"MAR02_SOURCE_UNIT_SUMMARY.tsv")
lab={z["unit_id"]:z["biological_label"] for z in labels}
_,allr=read_tsv(OUT/"MAR03_gene_set_results_all.tsv")
src=[]
for z in allr:
    zz=dict(z); zz["biological_label"]=lab.get(z["target_id"],z["target_id"]); src.append(zz)
write_tsv(OUT/"source_data"/"MAR03_SOURCE_GENE_SET_RESULTS.tsv",
          ["dataset","family","target_id","biological_label","model","frozen_n","mapped_n","coverage",
           "NGENES","beta","BETA_STD","SE","P_one","q_BH","primary_status"],src)

for f,outname in [
    ("MAR03_matched_meanZ.tsv","MAR03_SOURCE_MATCHED_CONTINUOUS.tsv"),
    ("MAR03_sensitivity.tsv","MAR03_SOURCE_SENSITIVITY.tsv"),
    ("MAR03_preanalysis_coverage.tsv","MAR03_SOURCE_COVERAGE.tsv"),
    ("MAR03_matching_summary.tsv","MAR03_SOURCE_MATCHING_SUMMARY.tsv"),
]:
    shutil.copy2(OUT/f,OUT/"source_data"/outname)

# Primary classification using only A marginal as headline; matched is complementary.
_,matched=read_tsv(OUT/"MAR03_matched_meanZ.tsv")
rows=[]
for t in sorted(lab):
    pg=next(z for z in allr if z["dataset"]=="PGC" and z["family"]=="A" and z["model"]=="MARGINAL" and z["target_id"]==t)
    sp=next(z for z in allr if z["dataset"]=="SPARK" and z["family"]=="A" and z["model"]=="MARGINAL" and z["target_id"]==t)
    mp=next(z for z in matched if z["dataset"]=="PGC" and z["family"]=="A" and z["target_id"]==t)
    ms=next(z for z in matched if z["dataset"]=="SPARK" and z["family"]=="A" and z["target_id"]==t)
    pgpos=float(pg["beta"])>0 and float(pg["q_BH"])<0.05
    sppos=float(sp["beta"])>0 and float(sp["q_BH"])<0.05
    matchpos=float(mp["q_BH_meanZ"])<0.05 and float(ms["q_BH_meanZ"])<0.05
    rows.append({
        "target_id":t,"biological_label":lab[t],
        "PGC_marginal_class":"POSITIVE" if pgpos else "NOT_DETECTED",
        "SPARK_marginal_class":"POSITIVE" if sppos else "NOT_DETECTED",
        "matched_meanZ_class":"BOTH_EMPIRICAL_POSITIVE" if matchpos else "NO_REPRODUCIBLE_EMPIRICAL_SIGNAL",
        "primary_reproducible_competitive_enrichment":"YES" if pgpos and sppos else "NO",
        "allowed_language":"reproducible common-variant enrichment detected" if pgpos and sppos
                           else "no reproducible common-variant enrichment detected among the prespecified higher-order units at the tested resolution"
    })
write_tsv(OUT/"MAR03_primary_classification.tsv",list(rows[0].keys()),rows)
shutil.copy2(OUT/"MAR03_primary_classification.tsv",OUT/"source_data"/"MAR03_SOURCE_PRIMARY_CLASSIFICATION.tsv")

# Locks/checksums/pass
result_files=[p for p in OUT.rglob("*") if p.is_file() and p.name not in {"MAR03_checksums.sha256","MAR03_PASS.txt"}]
lock={"stage":"MAR03","status":"PASS","version":"1.0R2_REAL","timestamp_utc":utc(),
      "preanalysis_lock_sha256":sha256_file(OUT/"MAR03_preanalysis_lock.json"),
      "output_hashes":{str(p.relative_to(OUT)):sha256_file(p) for p in result_files}}
(OUT/"MAR03_analysis_lock.json").write_text(json.dumps(lock,indent=2),encoding="utf-8")
lines=[f"{sha256_file(p)}  {p.relative_to(MA)}" for p in sorted([p for p in OUT.rglob("*") if p.is_file() and p.name!="MAR03_checksums.sha256"])]
(OUT/"MAR03_checksums.sha256").write_text("\n".join(lines)+"\n",encoding="utf-8")
(OUT/"MAR03_validation_summary.txt").write_text(
"""MAR03 technical execution PASS
Primary: 7 raw26-derived higher-order Expanded units, PGC/SPARK marginal competitive MAGMA.
Secondary completeness: 13 standalones.
Conditional: SUPPORTIVE_PENDING_MAR05 collinearity gate.
Matched continuous: historical MDV8 Match-A1 method; 10,000 unique balanced sets per A/B target; same sets scored in PGC and SPARK.
Sensitivities: Compact w10, Expanded remove-raw26 w10, Expanded w0/w50/MHC-excluded.
""",encoding="utf-8")
(OUT/"MAR03_PASS.txt").write_text(
"status=PASS\nstage=MAR03_COMMON_VARIANT_RETEST\nversion=1.0R2_REAL\n"
"primary_higher_order_n=7\nsecondary_standalone_n=13\nmatched_targets_n=20\nmatched_B=10000\n"
"conditional_status=SUPPORTIVE_PENDING_MAR05\n",encoding="utf-8")
print("MAR03_POSTFLIGHT=PASS")
