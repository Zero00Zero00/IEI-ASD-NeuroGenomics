
from common import *

ensure_dirs()
lockp=OUT/"MAR03_preanalysis_lock.json"; tokenp=OUT/"MAR03_STOPA_RELEASE.txt"
if not lockp.is_file() or not tokenp.is_file(): raise SystemExit("HOLD: freeze/release required")
tok={}
for line in tokenp.read_text().splitlines():
    if "=" in line:
        k,v=line.split("=",1); tok[k]=v
if tok.get("preanalysis_lock_sha256")!=sha256_file(lockp):
    raise SystemExit("HOLD: release token mismatch")

resolved=json.loads((OUT/"MAR03_resolved_inputs.json").read_text())
magma=Path(resolved["MAGMA_BIN"])
_,cov=read_tsv(OUT/"MAR03_preanalysis_coverage.tsv")

def run_one(ds,genes_raw,setfile,stype):
    targets=sorted({z["target_id"] for z in cov if z["dataset"]==ds and z["family"]=="A"})
    sdir=OUT/"sensitivity"/ds/stype
    sdir.mkdir(parents=True,exist_ok=True)
    pref=sdir/"marginal"
    cmd=[magma,"--gene-results",genes_raw,"--set-annot",setfile,
         "--model","direction=greater","--out",pref]
    run_cmd(cmd,str(pref)+".wrapper.log")
    _,rows=parse_gsa(gsa_output(pref))
    out=[]
    for t in targets:
        r=get_result_row(rows,t)
        ng,beta,bstd,se,p=result_values(r)
        out.append({"dataset":ds,"target_id":t,"sensitivity_type":stype,
                    "NGENES":ng,"beta":beta,"BETA_STD":bstd,"SE":se,"P_one":p,"q_BH":""})
    qs=bh([z["P_one"] for z in out])
    for z,q in zip(out,qs): z["q_BH"]=q
    return out

allr=[]
for ds,root,prefix in [("PGC",PGC,"PGC2019"),("SPARK",SPARK,"SPARK_EUR")]:
    allr += run_one(ds,Path(resolved[f"{ds}_w10_raw"]),OUT/"work"/f"MAR03_{ds}_compact.sets","COMPACT_W10")
    allr += run_one(ds,Path(resolved[f"{ds}_w10_raw"]),OUT/"work"/f"MAR03_{ds}_expanded_remove_raw26.sets","EXPANDED_REMOVE_RAW26_W10")
    allr += run_one(ds,Path(resolved[f"{ds}_w0_raw"]),OUT/"work"/f"MAR03_{ds}_higher_order_expanded.sets","EXPANDED_W0")
    allr += run_one(ds,Path(resolved[f"{ds}_w50_raw"]),OUT/"work"/f"MAR03_{ds}_higher_order_expanded.sets","EXPANDED_W50")
    allr += run_one(ds,Path(resolved[f"{ds}_mhc_raw"]),OUT/"work"/f"MAR03_{ds}_higher_order_expanded.sets","EXPANDED_MHCEXCLUDED")

write_tsv(OUT/"MAR03_sensitivity.tsv",
          ["dataset","target_id","sensitivity_type","NGENES","beta","BETA_STD","SE","P_one","q_BH"],allr)
print("MAR03_SENSITIVITY=PASS")
