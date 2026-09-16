
from common import *

ensure_dirs()
lockp=OUT/"MAR03_preanalysis_lock.json"; tokenp=OUT/"MAR03_STOPA_RELEASE.txt"
if not lockp.is_file() or not tokenp.is_file():
    raise SystemExit("HOLD: freeze + manual release required")
tok={}
for line in tokenp.read_text().splitlines():
    if "=" in line:
        k,v=line.split("=",1); tok[k]=v
if tok.get("preanalysis_lock_sha256")!=sha256_file(lockp):
    raise SystemExit("HOLD: release token does not match current preanalysis lock")
resolved=json.loads((OUT/"MAR03_resolved_inputs.json").read_text())
magma=Path(resolved["MAGMA_BIN"])

_,cov=read_tsv(OUT/"MAR03_preanalysis_coverage.tsv")
covmap={(z["dataset"],z["family"],z["target_id"]):z for z in cov}

def run_family(ds,genes_raw,setfile,fam,status,conditional=False):
    targets=sorted({z["target_id"] for z in cov if z["dataset"]==ds and z["family"]==fam})
    sdir=OUT/"magma"/ds/fam
    sdir.mkdir(parents=True,exist_ok=True)

    # Exact historical marginal CLI.
    pref=sdir/"marginal"
    cmd=[magma,"--gene-results",genes_raw,"--set-annot",setfile,
         "--model","direction=greater","--out",pref]
    run_cmd(cmd,str(pref)+".wrapper.log")
    _,rows=parse_gsa(gsa_output(pref))
    out=[]
    for t in targets:
        r=get_result_row(rows,t)
        ng,beta,bstd,se,p=result_values(r)
        c=covmap[(ds,fam,t)]
        out.append({
            "dataset":ds,"family":fam,"target_id":t,"model":"MARGINAL",
            "frozen_n":c["frozen_n"],"mapped_n":c["mapped_n"],"coverage":c["coverage"],
            "NGENES":ng,"beta":beta,"BETA_STD":bstd,"SE":se,"P_one":p,"q_BH":"",
            "primary_status":status
        })
    qs=bh([z["P_one"] for z in out])
    for z,q in zip(out,qs): z["q_BH"]=q

    if conditional and len(targets)>1:
        crows=[]
        for t in targets:
            others=[x for x in targets if x!=t]
            tf=sdir/f"target_{t}.txt"; tf.write_text(t+"\n",encoding="utf-8")
            pref=sdir/f"conditional_{t}"
            # Exact historical conditional CLI.
            cmd=[magma,"--gene-results",genes_raw,"--set-annot",setfile,
                 "--model",f"analyse=file,{tf}",f"condition-hide={','.join(others)}","direction=greater",
                 "--out",pref]
            run_cmd(cmd,str(pref)+".wrapper.log")
            _,rows=parse_gsa(gsa_output(pref))
            r=get_result_row(rows,t)
            ng,beta,bstd,se,p=result_values(r)
            c=covmap[(ds,fam,t)]
            crows.append({
                "dataset":ds,"family":fam,"target_id":t,"model":"CONDITIONAL_ALL_OTHER_A_UNITS",
                "frozen_n":c["frozen_n"],"mapped_n":c["mapped_n"],"coverage":c["coverage"],
                "NGENES":ng,"beta":beta,"BETA_STD":bstd,"SE":se,"P_one":p,"q_BH":"",
                "primary_status":"SUPPORTIVE_PENDING_MAR05"
            })
        qs=bh([z["P_one"] for z in crows])
        for z,q in zip(crows,qs): z["q_BH"]=q
        out += crows
    return out

datasets={
 "PGC":Path(resolved["PGC_w10_raw"]),
 "SPARK":Path(resolved["SPARK_w10_raw"]),
}
allr=[]
for ds,graw in datasets.items():
    allr += run_family(ds,graw,OUT/"work"/f"MAR03_{ds}_higher_order_expanded.sets","A","PRIMARY",True)
    allr += run_family(ds,graw,OUT/"work"/f"MAR03_{ds}_standalones.sets","B","SECONDARY_COMPLETENESS",False)
    allr += run_family(ds,graw,OUT/"work"/f"MAR03_{ds}_compact.sets","C","SENSITIVITY",False)

fields=["dataset","family","target_id","model","frozen_n","mapped_n","coverage","NGENES","beta","BETA_STD","SE","P_one","q_BH","primary_status"]
write_tsv(OUT/"MAR03_gene_set_results_all.tsv",fields,allr)

for ds in ["PGC","SPARK"]:
    rows=[z for z in allr if z["dataset"]==ds and z["family"]=="A" and z["model"]=="MARGINAL"]
    write_tsv(OUT/f"MAR03_{ds.lower()}_higher_order.tsv",fields,rows)
stand=[z for z in allr if z["family"]=="B" and z["model"]=="MARGINAL"]
write_tsv(OUT/"MAR03_standalone_tests.tsv",fields,stand)

(OUT/"state"/"MAR03_GSA_PASS.txt").write_text(
    "status=PASS\nmethod=historical_MDV6_MDV7_MAGMA_GSA\nmodel_direction=greater\n"
    "family_A_m=7\nfamily_B_m=13\nconditional=SUPPORTIVE_PENDING_MAR05\n",encoding="utf-8")
print("MAR03_GENE_SET_ANALYSIS=PASS")
