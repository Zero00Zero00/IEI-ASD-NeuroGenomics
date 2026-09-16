
from common import *

ensure_dirs()
checks=[]
def ck(cid,ok,obs,exp="",detail=""):
    checks.append({"check_id":cid,"status":"PASS" if ok else "HOLD","observed":obs,"expected":exp,"detail":detail})

req={
 "MAR02_PASS":MAR02/"MAR02_PASS.txt",
 "MAR02_COMMON_LOCK":MAR02/"MAR02_COMMON_HYPOTHESIS_LOCK.json",
 "MAR02_MEMBERSHIP":MAR02/"MAR02_raw26_domain_membership.tsv",
 "MAR02_STANDALONES":MAR02/"MAR02_raw26_standalones.tsv",
}
for k,p in req.items(): ck(k,p.is_file(),str(p),"exists")
if all(p.is_file() for p in req.values()):
    ck("MAR02_PASS_status","status=PASS" in req["MAR02_PASS"].read_text(),"status=PASS","status=PASS")
    lock=json.loads(req["MAR02_COMMON_LOCK"].read_text())
    ck("MAR02_no_GWAS_leak",lock.get("association_results_read")=="NO" and lock.get("gwas_inputs_used")=="NO",
       f"{lock.get('association_results_read')}/{lock.get('gwas_inputs_used')}","NO/NO")
    ck("MAR02_unit_n",len(lock["families"]["A_primary_higher_order_Expanded"])==7,
       len(lock["families"]["A_primary_higher_order_Expanded"]),7)
    ck("MAR02_standalone_n",len(lock["families"]["B_secondary_standalones"])==13,
       len(lock["families"]["B_secondary_standalones"]),13)

paths={
 "PGC_gene_results":PGC/"PGC2019_gene_results_w10.tsv",
 "PGC_gene_map":PGC/"MDV6_gene_id_map.tsv",
 "PGC_PASS":PGC/"MDV6_PASS.txt",
 "PGC_w10_raw":PGC/"magma"/"PGC2019_w10.genes.raw",
 "PGC_w0_raw":PGC/"magma"/"PGC2019_w0.genes.raw",
 "PGC_w50_raw":PGC/"magma"/"PGC2019_w50.genes.raw",
 "PGC_mhc_raw":PGC/"magma"/"PGC2019_mhc.genes.raw",
 "SPARK_gene_results":SPARK/"SPARK_EUR_gene_results_w10.tsv",
 "SPARK_gene_map":SPARK/"MDV7_gene_id_map.tsv",
 "SPARK_PASS":SPARK/"MDV7_PASS.txt",
 "SPARK_w10_raw":SPARK/"magma"/"SPARK_EUR_w10.genes.raw",
 "SPARK_w0_raw":SPARK/"magma"/"SPARK_EUR_w0.genes.raw",
 "SPARK_w50_raw":SPARK/"magma"/"SPARK_EUR_w50.genes.raw",
 "SPARK_mhc_raw":SPARK/"magma"/"SPARK_EUR_mhc.genes.raw",
 "MDV8_common":MDV8/"MDV8_common_scored_gene_universe.tsv",
 "MDV8_PASS":MDV8/"MDV8_PASS.txt",
 "MDV6_runner":CH3/"run_mdv6_staged.sh",
 "MDV7_runner":CH3/"run_mdv7_staged.sh",
 "MDV6_gsa_script":CH3/"workflow"/"scripts"/"mdv6_05_run_gsa.py",
 "MDV7_gsa_script":CH3/"workflow"/"scripts"/"mdv7_05_run_gsa.py",
 "MDV8_match_script":CH3/"workflow"/"scripts"/"mdv8_05_matched_null.py",
 "MDV8_config":CH3/"config"/"mdv8_specificity_gate_v1_1.yaml",
}
for k,p in paths.items(): ck(k,p.is_file(),str(p),"exists")

# frozen authority checksum checks for main w10 inputs
for label,target,manifest in [
 ("PGC_gene_results_checksum",paths["PGC_gene_results"],PGC/"MDV6_checksums.sha256"),
 ("PGC_w10_checksum",paths["PGC_w10_raw"],PGC/"MDV6_checksums.sha256"),
 ("PGC_w0_checksum",paths["PGC_w0_raw"],PGC/"MDV6_checksums.sha256"),
 ("PGC_w50_checksum",paths["PGC_w50_raw"],PGC/"MDV6_checksums.sha256"),
 ("PGC_mhc_checksum",paths["PGC_mhc_raw"],PGC/"MDV6_checksums.sha256"),
 ("SPARK_gene_results_checksum",paths["SPARK_gene_results"],SPARK/"MDV7_checksums.sha256"),
 ("SPARK_w10_checksum",paths["SPARK_w10_raw"],SPARK/"MDV7_checksums.sha256"),
 ("SPARK_w0_checksum",paths["SPARK_w0_raw"],SPARK/"MDV7_checksums.sha256"),
 ("SPARK_w50_checksum",paths["SPARK_w50_raw"],SPARK/"MDV7_checksums.sha256"),
 ("SPARK_mhc_checksum",paths["SPARK_mhc_raw"],SPARK/"MDV7_checksums.sha256"),
]:
    st,exp=verify_manifest_entry(manifest,target); ck(label,st=="PASS",st,"PASS",f"expected={exp}")

# Real historical interface semantics
if paths["MDV6_gsa_script"].is_file():
    t=paths["MDV6_gsa_script"].read_text(encoding="utf-8",errors="replace")
    ok="--model" in t and "direction=" in t and "condition-hide=" in t and "analyse=file" in t
    ck("MDV6_GSA_interface_semantics",ok,"historical --model/direction/conditional syntax","present")
if paths["MDV7_gsa_script"].is_file():
    t=paths["MDV7_gsa_script"].read_text(encoding="utf-8",errors="replace")
    ok="--model" in t and "direction=" in t and "condition-hide=" in t and "analyse=file" in t
    ck("MDV7_GSA_interface_semantics",ok,"historical --model/direction/conditional syntax","present")
if paths["MDV8_match_script"].is_file() and paths["MDV8_config"].is_file():
    st=paths["MDV8_match_script"].read_text(encoding="utf-8",errors="replace")
    ct=paths["MDV8_config"].read_text(encoding="utf-8",errors="replace")
    ok=("chromosome_stratified_exponential_tilt_rerandomization" in st and
        "exact_chromosome_composition" in ct and
        "log_gene_length_z" in ct and "GC_z_mdv8" in ct and "log1p_annotation_degree_z" in ct)
    ck("MDV8_matching_interface_semantics",ok,"historical exponential-tilt exact-chromosome matching","present")

try:
    magma,allm=find_magma()
    cp=subprocess.run([str(magma),"--version"],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=30)
    ck("MAGMA",cp.returncode==0,str(magma),"executable",cp.stdout.strip()[:300])
except Exception as e:
    magma=None; ck("MAGMA",False,repr(e),"executable")

# validated Python direct runtime
if PYTHON_BIN.is_file():
    cp=subprocess.run([str(PYTHON_BIN),"-c","import numpy,scipy,yaml; print('OK')"],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    ck("PYTHON_RUNTIME",cp.returncode==0,str(PYTHON_BIN),"validated executable",cp.stdout.strip())
else:
    ck("PYTHON_RUNTIME",False,str(PYTHON_BIN),"exists")

snap=[]
for k in ["MDV6_runner","MDV7_runner","MDV6_gsa_script","MDV7_gsa_script","MDV8_match_script","MDV8_config"]:
    p=paths[k]; snap.append({"role":k,"path":str(p),"sha256":sha256_file(p) if p.is_file() else "MISSING"})
write_tsv(OUT/"MAR03_HISTORICAL_INTERFACE_SNAPSHOT.tsv",["role","path","sha256"],snap)

resolved={k:str(v) for k,v in paths.items()}
resolved.update({"MAGMA_BIN":str(magma) if magma else "","PYTHON_BIN":str(PYTHON_BIN)})
(OUT/"MAR03_resolved_inputs.json").write_text(json.dumps(resolved,indent=2),encoding="utf-8")
write_tsv(OUT/"MAR03_preflight_report.tsv",["check_id","status","observed","expected","detail"],checks)
bad=[z for z in checks if z["status"]!="PASS"]
if bad:
    print(f"HOLD: MAR03 R2 preflight failed ({len(bad)} checks). See {OUT/'MAR03_preflight_report.tsv'}")
    raise SystemExit(20)
print("MAR03_PREFLIGHT=PASS")
print("METHOD_AUTHORITY=MDV6/MDV7_GSA + MDV8_MATCH_A1")
