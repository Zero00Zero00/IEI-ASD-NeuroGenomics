
from common import *

ensure_dirs()
pre=OUT/"MAR03_preflight_report.tsv"
if not pre.is_file(): raise SystemExit("HOLD: run preflight first")
_,pre_rows=read_tsv(pre)
if any(z["status"]!="PASS" for z in pre_rows): raise SystemExit("HOLD: preflight not PASS")
resolved=json.loads((OUT/"MAR03_resolved_inputs.json").read_text())

mh,mrows=read_tsv(MAR02/"MAR02_raw26_domain_membership.tsv")
sh,srows=read_tsv(MAR02/"MAR02_raw26_standalones.tsv")
_,labelrows=read_tsv(PKG/"config"/"MAR02_SOURCE_UNIT_SUMMARY.tsv")
labels={z["unit_id"]:z["biological_label"] for z in labelrows}

# A and C membership from frozen MAR02.
members={"A":{},"C":{},"B":{}}
for z in mrows:
    u=z["unit_id"]; sym=(z.get("gene") or "").strip().upper(); hg=(z.get("HGNC_id") or "").strip().upper()
    ident=(sym,hg)
    if z.get("Expanded_flag")=="1": members["A"].setdefault(u,set()).add(ident)
    if z.get("Compact_flag")=="1": members["C"].setdefault(u,set()).add(ident)

# B membership from frozen standalone pathway keys + global membership.
pm_path,ph,prows,source_col,pid_col,gene_col,hgnc_col,key_col=resolve_pathway_membership()
stand_keys={z["standalone_id"]:z["pathway_key"] for z in srows}
for sid in stand_keys: members["B"][sid]=set()
for sid,pkey in stand_keys.items():
    src,pid=pkey.split("::",1)
    for z in prows:
        match=False
        if key_col and (z.get(key_col) or "").strip()==pkey: match=True
        elif source_col and pid_col and (z.get(source_col) or "").strip()==src and (z.get(pid_col) or "").strip()==pid: match=True
        if not match: continue
        sym=(z.get(gene_col) or "").strip().upper() if gene_col else ""
        hg=(z.get(hgnc_col) or "").strip().upper() if hgnc_col else ""
        if sym or hg: members["B"][sid].add((sym,hg))

counts={k:len(v) for k,v in members.items()}
if counts!={"A":7,"B":13,"C":7}:
    raise SystemExit(f"HOLD: family counts A/B/C={counts['A']}/{counts['B']}/{counts['C']}")
empty=[(fam,t) for fam in members for t,s in members[fam].items() if not s]
if empty: raise SystemExit("HOLD: empty frozen membership targets: "+str(empty))

maps={"PGC":detect_map(PGC/"MDV6_gene_id_map.tsv"),"SPARK":detect_map(SPARK/"MDV7_gene_id_map.tsv")}
membership_rows=[]; coverage_rows=[]

def map_members(ds,fam,target):
    m=maps[ds]; mapped=[]; missing=[]
    for sym,hg in sorted(members[fam][target]):
        gid=m["by_symbol"].get(sym) if sym else None
        if not gid and hg: gid=m["by_hgnc"].get(hg)
        if gid: mapped.append(gid)
        else: missing.append(sym or hg)
    return sorted(set(mapped),key=lambda x:(not str(x).isdigit(),str(x))),sorted(set(missing))

for ds in ["PGC","SPARK"]:
    for fam in ["A","B","C"]:
        for target in sorted(members[fam]):
            gids,missing=map_members(ds,fam,target)
            frozen_n=len(members[fam][target])
            coverage_rows.append({
                "dataset":ds,"family":fam,"target_id":target,"frozen_n":frozen_n,
                "mapped_n":len(gids),"coverage":len(gids)/frozen_n if frozen_n else "",
                "missing_n":len(missing),"missing_genes":";".join(missing)
            })
            for gid in gids:
                membership_rows.append({"dataset":ds,"family":fam,"target_id":target,"magma_gene_id":gid})

write_tsv(OUT/"MAR03_frozen_membership.tsv",["dataset","family","target_id","magma_gene_id"],membership_rows)
write_tsv(OUT/"MAR03_preanalysis_coverage.tsv",
          ["dataset","family","target_id","frozen_n","mapped_n","coverage","missing_n","missing_genes"],coverage_rows)

# Strict identity check: same gene-id maps are expected and observed in the real probe.
for fam in ["A","B","C"]:
    for target in members[fam]:
        pg={z["magma_gene_id"] for z in membership_rows if z["dataset"]=="PGC" and z["family"]==fam and z["target_id"]==target}
        sp={z["magma_gene_id"] for z in membership_rows if z["dataset"]=="SPARK" and z["family"]==fam and z["target_id"]==target}
        if pg!=sp:
            raise SystemExit(f"HOLD: PGC/SPARK target gene IDs differ {fam}/{target}")

# Build MAGMA .sets files.
for ds in ["PGC","SPARK"]:
    for fam,name in [("A","higher_order_expanded"),("B","standalones"),("C","compact")]:
        p=OUT/"work"/f"MAR03_{ds}_{name}.sets"
        with open(p,"w",encoding="utf-8",newline="\n") as f:
            for target in sorted(members[fam]):
                gids=[z["magma_gene_id"] for z in membership_rows if z["dataset"]==ds and z["family"]==fam and z["target_id"]==target]
                f.write(target+" "+" ".join(gids)+"\n")

# remove-raw26 A sensitivity
raw26_syms={(z.get("gene") or "").strip().upper() for z in mrows if z.get("raw26_anchor_flag")=="1"}
for ds in ["PGC","SPARK"]:
    mp=maps[ds]
    anchor_ids={mp["by_symbol"].get(s) for s in raw26_syms if mp["by_symbol"].get(s)}
    p=OUT/"work"/f"MAR03_{ds}_expanded_remove_raw26.sets"
    with open(p,"w",encoding="utf-8",newline="\n") as f:
        for target in sorted(members["A"]):
            gids=[z["magma_gene_id"] for z in membership_rows if z["dataset"]==ds and z["family"]=="A" and z["target_id"]==target]
            gids=sorted(set(gids)-set(anchor_ids),key=lambda x:(not str(x).isdigit(),str(x)))
            f.write(target+" "+" ".join(gids)+"\n")

# canonical hypothesis snapshot
_,hrows=read_tsv(PKG/"config"/"MAR02_SOURCE_MAR03_HYPOTHESES.tsv")
write_tsv(OUT/"MAR03_frozen_hypotheses.tsv",list(hrows[0].keys()),hrows)

# Freeze before any association-score parsing.
authority=[
 Path(resolved["PGC_gene_results"]),Path(resolved["PGC_w10_raw"]),Path(resolved["PGC_gene_map"]),
 Path(resolved["SPARK_gene_results"]),Path(resolved["SPARK_w10_raw"]),Path(resolved["SPARK_gene_map"]),
 Path(resolved["MDV8_common"]),MAR02/"MAR02_COMMON_HYPOTHESIS_LOCK.json",
 MAR02/"MAR02_raw26_domain_membership.tsv",MAR02/"MAR02_raw26_standalones.tsv",
 Path(resolved["MDV6_gsa_script"]),Path(resolved["MDV7_gsa_script"]),
 Path(resolved["MDV8_match_script"]),Path(resolved["MDV8_config"])
]
lock={
 "stage":"MAR03_COMMON_VARIANT_RETEST",
 "version":"1.0R2_REAL",
 "timestamp_utc":utc(),
 "association_scores_parsed":"NO",
 "magma_gene_window_kb":10,
 "families":{"A_primary_higher_order_expanded_n":7,"B_secondary_standalone_n":13,"C_compact_sensitivity_n":7},
 "primary_model":"MAGMA competitive marginal; historical MDV6/MDV7 CLI with --model direction=greater",
 "conditional_model":"historical per-target analyse=file + condition-hide all other A units; SUPPORTIVE_PENDING_MAR05",
 "matched_continuous":{
   "targets":"A+B",
   "meanZ":"PRIMARY_COMPLEMENTARY",
   "medianZ":"SUPPORTIVE",
   "same_matched_sets_PGC_SPARK":"YES",
   "B":10000,
   "seed":20260814,
   "matching_method":"chromosome_stratified_exponential_tilt_rerandomization",
   "exact_chromosome_composition":"YES",
   "distance_variables":["log_gene_length_z","GC_z_mdv8","log1p_annotation_degree_z"],
   "sampling_weight":"exponential_tilt_on_frozen_structural_covariates",
   "per_set_abs_SMD_caliper":0.25,
   "median_abs_SMD_max":0.15,
   "p95_abs_SMD_max":0.25,
   "tilt_mean_residual_max":0.02,
   "max_set_generation_attempt_factor":40,
   "control_exclusion":"ALL_FROZEN_MAR03_A_PLUS_B_HYPOTHESIS_GENES",
   "historical_method_authority":"MDV8_v1_1 mdv8_05_matched_null.py + mdv8_specificity_gate_v1_1.yaml",
   "extension_note":"Historical MDV8 matching algorithm is prospectively extended to the secondary standalone family before any association score was read."
 },
 "sensitivity":["Compact_w10","Expanded_remove_raw26_w10","Expanded_w0","Expanded_w50","Expanded_MHCexcluded"],
 "multiplicity":{"A_m":7,"B_m":13,"C_m":7,"matched_BH":"separate within dataset × family × statistic"},
 "input_hashes":{str(p):sha256_file(p) for p in authority},
 "set_file_hashes":{p.name:sha256_file(p) for p in sorted((OUT/"work").glob("*.sets"))},
 "GWAS_score_selection_used":"NO"
}
(OUT/"MAR03_preanalysis_lock.json").write_text(json.dumps(lock,indent=2),encoding="utf-8")

# Explicit method amendment record (R1 was never released).
(OUT/"MAR03_METHOD_IMPLEMENTATION_NOTE.txt").write_text(
"""status=FROZEN_PRE_ASSOCIATION
reason=Real historical interface probe showed that the R1 provisional matching contract did not match frozen MDV8 Match-A1.
action=Use historical chromosome-stratified exponential-tilt rerandomization with exact chromosome composition and MDV8 frozen structural covariates.
scientific_hypotheses_changed=NO
association_scores_read=NO
family_A=7
family_B=13
family_C=7
""",encoding="utf-8")

print("MAR03_FREEZE=PASS")
print("A_B_C=7/13/7")
print("MANUAL_GATE_REQUIRED=YES")
