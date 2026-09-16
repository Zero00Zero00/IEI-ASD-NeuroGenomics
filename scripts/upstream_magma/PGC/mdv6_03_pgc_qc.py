#!/usr/bin/env python3
import argparse, gzip, math, os, json
from collections import Counter
from pathlib import Path
import pandas as pd
from mdv6_common import *

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); a=ap.parse_args()
    cfg=load_cfg(a.config); root=Path(cfg["project_root"]); out=out_dir(cfg)
    require_file(out/cfg["outputs"]["preanalysis_pass"])
    resolved=json.load(open(out/cfg["outputs"]["resolved_inputs"]))
    pgc=Path(resolved["pgc_gwas"]); bim=resolved["ld_prefix"]+".bim"
    # reference dictionary
    ref={}; refdup=set(); ref_rows=0
    with open(bim,encoding="utf-8",errors="replace") as f:
        for line in f:
            a0=line.split();
            if len(a0)<6: continue
            ch=norm_chr(a0[0]); snp=a0[1]; bp=int(float(a0[3])); a1=a0[4]; a2=a0[5]; ref_rows+=1
            if snp in ref: refdup.add(snp)
            else: ref[snp]=(ch,bp,a1,a2)
    for x in refdup: ref.pop(x,None)
    # count PGC duplicate IDs in a first pass using pandas chunks
    counts=Counter(); total=0
    for ch in pd.read_csv(pgc,sep=r"\s+",usecols=["SNP"],dtype={"SNP":str},chunksize=700000,compression="gzip"):
        counts.update(ch["SNP"].astype(str)); total+=len(ch)
    dupids={k for k,v in counts.items() if v>1}
    # second pass, QC and stream outputs
    qcgz=out/cfg["outputs"]["pgc_qc"]; pval=out/cfg["outputs"]["pgc_pval"]; sloc=out/cfg["outputs"]["pgc_snploc"]
    pvalm=out/cfg["outputs"]["pgc_pval_mhc"]; slocm=out/cfg["outputs"]["pgc_snploc_mhc"]
    reason=Counter(); allele_status=Counter(); kept=0; kept_mhc=0
    with gzip.open(qcgz,"wt",encoding="utf-8") as fq, open(pval,"w") as fp, open(sloc,"w") as fs, open(pvalm,"w") as fpm, open(slocm,"w") as fsm:
        fq.write("SNP\tCHR\tBP\tA1\tA2\tINFO\tP\tREF_A1\tREF_A2\tallele_status\n"); fp.write("SNP\tP\n"); fs.write("SNP\tCHR\tBP\n"); fpm.write("SNP\tP\n"); fsm.write("SNP\tCHR\tBP\n")
        use=["CHR","SNP","BP","A1","A2","INFO","P"]
        for df in pd.read_csv(pgc,sep=r"\s+",usecols=use,chunksize=300000,compression="gzip",dtype={"SNP":str,"A1":str,"A2":str}):
            for r in df.itertuples(index=False):
                snp=str(r.SNP); ch=norm_chr(r.CHR)
                try: bp=int(float(r.BP)); info=float(r.INFO); pv=float(r.P)
                except: reason["NONNUMERIC"]+=1; continue
                if snp in dupids: reason["PGC_DUPLICATE_SNP_ID"]+=1; continue
                if not (math.isfinite(pv) and 0<pv<=1): reason["INVALID_P"]+=1; continue
                if not (math.isfinite(info) and info>=float(cfg["pgc2019"]["info_min_primary"])): reason["INFO_BELOW_THRESHOLD"]+=1; continue
                if ch is None or bp<=0: reason["INVALID_POSITION"]+=1; continue
                rr=ref.get(snp)
                if rr is None: reason["NOT_UNIQUE_IN_1KG_REFERENCE"]+=1; continue
                rch,rbp,ra1,ra2=rr
                if cfg["pgc2019"]["require_reference_chr_bp_match"] and (ch!=rch or bp!=rbp): reason["CHR_BP_MISMATCH"]+=1; continue
                ok,ast=allele_compatible(r.A1,r.A2,ra1,ra2); allele_status[ast]+=1
                if cfg["pgc2019"]["require_allele_compatible"] and not ok: reason["ALLELE_MISMATCH"]+=1; continue
                fq.write(f"{snp}\t{ch}\t{bp}\t{r.A1}\t{r.A2}\t{info:.6g}\t{pv:.12g}\t{ra1}\t{ra2}\t{ast}\n")
                fp.write(f"{snp}\t{pv:.12g}\n"); fs.write(f"{snp}\t{ch}\t{bp}\n"); kept+=1
                in_mhc=(ch==int(cfg["sensitivity"]["mhc_chr"]) and int(cfg["sensitivity"]["mhc_start_bp"])<=bp<=int(cfg["sensitivity"]["mhc_end_bp"]))
                if not in_mhc:
                    fpm.write(f"{snp}\t{pv:.12g}\n"); fsm.write(f"{snp}\t{ch}\t{bp}\n"); kept_mhc+=1
                else: reason["MHC_EXCLUDED_SENSITIVITY"]+=1
    report=[]
    def add(metric,value,expected="",status="PASS"): report.append({"metric":metric,"value":value,"expected_or_rule":expected,"status":status})
    add("PGC_input_rows",total,f">={cfg['pgc2019']['expected_min_rows']}","PASS" if total>=cfg["pgc2019"]["expected_min_rows"] else "FAIL")
    add("PGC_duplicate_unique_ids",len(dupids),cfg["pgc2019"]["duplicate_snp_rule"])
    add("1KG_reference_rows",ref_rows)
    add("1KG_nonunique_ids",len(refdup),"excluded")
    add("INFO_min",cfg["pgc2019"]["info_min_primary"])
    add("kept_reference_compatible_snps",kept,">1000000","PASS" if kept>1000000 else "FAIL")
    add("kept_MHCexcluded_snps",kept_mhc)
    add("sample_size_primary",cfg["pgc2019"]["sample_size_primary"],"18381+27969")
    for k,v in sorted(reason.items()): add("drop_"+k,v)
    for k,v in sorted(allele_status.items()): add("allele_"+k,v)
    rdf=pd.DataFrame(report); rdf.to_csv(out/cfg["outputs"]["pgc_qc_report"],sep="\t",index=False)
    if (rdf.status=="FAIL").any(): raise RuntimeError("PGC QC hard gate failed")
    write_text(out/cfg["outputs"]["pgc_qc_pass"],f"status=PASS\nstage=MDV6_PGC_QC\nassociation_results_read=YES\ninput_rows={total}\nkept_snps={kept}\nINFO_min={cfg['pgc2019']['info_min_primary']}\nreference_chr_bp_match=YES\nallele_compatible=YES\n")
    print("MDV6_PGC_QC_PASS",kept)
if __name__=="__main__": main()
