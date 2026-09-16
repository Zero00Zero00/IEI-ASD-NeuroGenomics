
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

ensure_dirs()
p=OUT/"MAR07_preflight_report.tsv"
if not p.is_file(): raise SystemExit("HOLD: preflight first")
_,rows=read_tsv(p)
if any(z["status"]!="PASS" for z in rows): raise SystemExit("HOLD: preflight not PASS")
resolved=json.loads((OUT/"MAR07_resolved_inputs.json").read_text())
lock={
 "stage":"MAR07_CHAIN_CALIBRATION_R4_CERTIFICATION",
 "version":"1.0R4.1",
 "timestamp_utc":utc(),
 "observed_authority":{
   "MDV2":"full frozen 6671 historical MDV2 results; 18 primary-pass crossstage rows",
   "MDV3":"WP01 raw26 exact legacy solver, 6671 pathways",
   "MDV4":"WP01 empirical tiering on frozen 10,000 matched sets",
   "Tier1":"MDV3 enriched AND q_Firth<0.05 AND q_emp<0.05; MDV2 not prerequisite"
 },
 "certification":{
   "MDV2":"primary-only adapter must numerically replay frozen full MDV2 primary results",
   "MDV3":"exact wp01_mdv3_legacy_bridge + mdv3_firth_solver_v1p3 replay",
   "MDV4":"exact wp01_empirical_tiering.empirical_for_sets replay",
   "numeric_tolerance_MDV2":1e-8,
   "numeric_tolerance_MDV3":1e-8,
   "numeric_tolerance_MDV4":1e-12,
   "required_replay_n":{"MDV2":6671,"MDV3":6671,"MDV4":6671},
   "numeric_missing_semantics":"both NA/NaN/blank = semantically equal; one-sided missing = FAIL; finite values use frozen tolerances"
 },
 "global_null":{
   "gene_universe_n":19267,
   "pathways_n":6671,
   "group_sizes":{"Overlap_raw":26,"IEI_only":474,"ASD_only":909,"Background":17858},
   "group_sampler":"joint-stratum exponential-tilt extension of frozen WP01 global matching engine",
   "joint_bins_per_covariate":5,
   "features":["z_log_gene_length","z_GC","z_log_annotation"],
   "per_group_abs_SMD_caliper":0.20,
   "group_order":["Overlap_raw","IEI_only","ASD_only"],
   "mutually_exclusive":"YES",
   "MDV2":"certified primary-only algorithm replay-equivalent to frozen MDV2 primary",
   "MDV3":"certified legacy solver",
   "MDV4":"frozen raw26 matched-set bank; remove any reference set intersecting pseudo-overlap genes",
   "pilot_R":50,"primary_R":250,"adaptive_R":[500,1000],
   "CI_halfwidth_target":0.02,
   "submission_HOLD_if_MC_CI_upper_gt":0.10,
   "seed_base":20260909
 },
 "no_rule_tuning_after_release":"YES",
 "input_hashes":resolved["hashes"]
}
(OUT/"MAR07_preanalysis_lock.json").write_text(json.dumps(lock,indent=2),encoding="utf-8")
print("MAR07_R4_FREEZE=PASS")
print("HOLD_FOR_MANUAL_GATE")
