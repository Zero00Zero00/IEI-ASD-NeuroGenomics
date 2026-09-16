#!/usr/bin/env python3
import argparse, gzip, math, json, os, re, subprocess, tempfile
from collections import Counter, defaultdict
from pathlib import Path
import pandas as pd
from mdv7_common import *

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); a=ap.parse_args(); cfg=load_cfg(a.config); root=Path(cfg["project_root"]); out=out_dir(cfg); prov=prov_dir(cfg)
    stopa=prov/"MDV7_STOPA_RELEASE.txt"; require_file(stopa,"STOP A release")
    resolved=json.load(open(out/cfg["outputs"]["resolved_inputs"])); sp=Path(resolved["spark_primary"]); liftover=resolved["liftover_bin"]; chain=resolved["liftover_chain"]; bim=resolved["ld_prefix"]+".bim"
    # first pass: ADD-row ID counts and total rows; association values are now legitimately read
    counts=Counter(); raw_rows=add_rows=0
    use=["CHROM","POS","ID","REF","ALT","A1","A2","TEST","OBS_CT","P"]
    for df in pd.read_csv(sp,sep=r"\s+",usecols=use,chunksize=400000,compression="gzip",dtype={"ID":str,"REF":str,"ALT":str,"A1":str,"A2":str,"TEST":str}):
        raw_rows+=len(df); ad=df[df.TEST.astype(str)==cfg["spark"]["additive_test_value"]]; add_rows+=len(ad); counts.update(ad.ID.astype(str))
    if raw_rows<int(cfg["spark"]["min_raw_rows"]): raise RuntimeError(f"SPARK file implausibly small: {raw_rows}")
    dup_raw={k for k,v in counts.items() if v>1}
    # second pass: pre-liftover technical QC and BED creation. Embed all required row info in BED name.
    tmp=out/"harmonization_tmp"; tmp.mkdir(exist_ok=True); bed=tmp/"spark_grch38.bed"; reasons=Counter(); id_checked=id_match=0; eligible=0
    with open(bed,"w",encoding="utf-8") as fb:
        for df in pd.read_csv(sp,sep=r"\s+",usecols=use,chunksize=250000,compression="gzip",dtype={"ID":str,"REF":str,"ALT":str,"A1":str,"A2":str,"TEST":str}):
            for r in df.itertuples(index=False):
                if str(r.TEST)!=cfg["spark"]["additive_test_value"]: reasons["NON_ADD_TEST"]+=1; continue
                vid=str(r.ID)
                if vid in dup_raw: reasons["DUPLICATE_RAW_ID"]+=1; continue
                ch=norm_chr(r.CHROM)
                try: bp=int(float(r.POS)); p=float(r.P); n=int(float(r.OBS_CT))
                except: reasons["NONNUMERIC"]+=1; continue
                if cfg["harmonization"]["autosomes_only"] and (ch is None or ch<1 or ch>22): reasons["NON_AUTOSOMAL"]+=1; continue
                if not (math.isfinite(p) and 0<p<=1): reasons["INVALID_P"]+=1; continue
                if n<=0: reasons["INVALID_OBS_CT"]+=1; continue
                emb=parse_spark_embedded_id(vid)
                if emb:
                    id_checked+=1; ech,ep,er,ea=emb
                    same=(ech==ch and ep==bp and er==str(r.REF).upper() and ea==str(r.ALT).upper())
                    if same: id_match+=1
                    else: reasons["EMBEDDED_ID_MISMATCH"]+=1; continue
                else: reasons["UNPARSEABLE_EMBEDDED_ID"]+=1; continue
                ok,_=allele_compatible(r.REF,r.ALT,r.A1,r.A2)
                if not ok: reasons["RAW_REFALT_A1A2_MISMATCH"]+=1; continue
                name="|".join([vid,str(r.REF),str(r.ALT),str(r.A1),str(r.A2),f"{p:.16g}",str(n),str(ch),str(bp)])
                fb.write(f"chr{ch}\t{bp-1}\t{bp}\t{name}\n"); eligible+=1
    id_rate=(id_match/id_checked) if id_checked else 0
    if id_rate<float(cfg["harmonization"]["embedded_id_match_min"]): raise RuntimeError(f"Embedded ID consistency {id_rate:.4%} below frozen gate")
    # liftOver all eligible GRCh38 positions to GRCh37
    lifted=tmp/"spark_grch37.bed"; unmapped=tmp/"spark_unmapped.bed"
    run([liftover,f"-minMatch={cfg['liftover']['min_match']}",bed,chain,lifted,unmapped],tmp/"liftOver.wrapper.log")
    lifted_n=sum(1 for _ in open(lifted,encoding="utf-8",errors="replace")); lift_rate=lifted_n/eligible if eligible else 0
    if lift_rate<float(cfg["harmonization"]["liftover_success_min"]): raise RuntimeError(f"liftOver success {lift_rate:.4%} below frozen gate")
    # build 1KG coordinate index and remove nonunique reference rsIDs
    idcounts=Counter(); ref_by_pos=defaultdict(list); ref_rows=0
    with open(bim,encoding="utf-8",errors="replace") as f:
        for line in f:
            z=line.split();
            if len(z)<6: continue
            ch=norm_chr(z[0]); sid=z[1]; bp=int(float(z[3])); a1=z[4]; a2=z[5]; ref_rows+=1; idcounts[sid]+=1; ref_by_pos[(ch,bp)].append((sid,a1,a2))
    # first candidate pass: lifted coordinate + unique compatible 1KG allele match
    cand=tmp/"candidate.tsv"; cand_counts=Counter(); candidate_n=0
    with open(lifted,encoding="utf-8",errors="replace") as fi, open(cand,"w",encoding="utf-8") as fo:
        fo.write("SNP\tP\tN\tCHR\tBP\tRAW_ID\tRAW_CHR\tRAW_BP\tREF\tALT\tA1\tA2\tREF_A1\tREF_A2\tALLELE_STATUS\n")
        for line in fi:
            z=line.rstrip("\n").split("\t");
            if len(z)<4: continue
            lch=norm_chr(z[0]); lbp=int(z[1])+1; parts=z[3].split("|")
            if len(parts)!=9: reasons["BAD_LIFTED_NAME"]+=1; continue
            vid,ref,alt,a1,a2,p,n,rawch,rawbp=parts
            hits=[]
            for sid,ra1,ra2 in ref_by_pos.get((lch,lbp),[]):
                if idcounts[sid]!=1: continue
                ok,ast=allele_compatible(ref,alt,ra1,ra2)
                if ok: hits.append((sid,ra1,ra2,ast))
            if len(hits)==0: reasons["NO_1KG_COORD_ALLELE_MATCH"]+=1; continue
            if len(hits)>1: reasons["AMBIGUOUS_1KG_COORD_ALLELE_MATCH"]+=1; continue
            sid,ra1,ra2,ast=hits[0]; candidate_n+=1; cand_counts[sid]+=1
            fo.write(f"{sid}\t{p}\t{n}\t{lch}\t{lbp}\t{vid}\t{rawch}\t{rawbp}\t{ref}\t{alt}\t{a1}\t{a2}\t{ra1}\t{ra2}\t{ast}\n")
    dup_refmap={k for k,v in cand_counts.items() if v>1}
    # final outputs with 1KG IDs + Build37 coords
    qcgz=out/cfg["outputs"]["harmonized_qc"]; pval=out/cfg["outputs"]["spark_pval"]; sloc=out/cfg["outputs"]["spark_snploc"]; pvalm=out/cfg["outputs"]["spark_pval_mhc"]; slocm=out/cfg["outputs"]["spark_snploc_mhc"]
    kept=kept_mhc=0; allele_status=Counter();
    with open(cand,encoding="utf-8") as fi, gzip.open(qcgz,"wt",encoding="utf-8") as fq, open(pval,"w") as fp, open(sloc,"w") as fs, open(pvalm,"w") as fpm, open(slocm,"w") as fsm:
        hdr=fi.readline(); fq.write(hdr); fp.write("SNP\tP\tN\n"); fs.write("SNP\tCHR\tBP\n"); fpm.write("SNP\tP\tN\n"); fsm.write("SNP\tCHR\tBP\n")
        for line in fi:
            z=line.rstrip("\n").split("\t"); sid=z[0]
            if sid in dup_refmap: reasons["DUPLICATE_FINAL_1KG_ID"]+=1; continue
            ch=int(z[3]); bp=int(z[4]); p=z[1]; n=z[2]; ast=z[-1]; allele_status[ast]+=1
            fq.write(line); fp.write(f"{sid}\t{p}\t{n}\n"); fs.write(f"{sid}\t{ch}\t{bp}\n"); kept+=1
            inmhc=(ch==int(cfg["sensitivity"]["mhc_chr"]) and int(cfg["sensitivity"]["mhc_start_bp"])<=bp<=int(cfg["sensitivity"]["mhc_end_bp"]))
            if not inmhc: fpm.write(f"{sid}\t{p}\t{n}\n"); fsm.write(f"{sid}\t{ch}\t{bp}\n"); kept_mhc+=1
    retention=kept/add_rows if add_rows else 0
    if retention<float(cfg["harmonization"]["reference_compatible_retention_min"]): raise RuntimeError(f"Reference-compatible retention {retention:.4%} below catastrophic-failure gate")
    # audit/report
    rows=[
      ("raw_rows",raw_rows),("additive_rows",add_rows),("duplicate_raw_ids",len(dup_raw)),("eligible_for_liftover",eligible),("embedded_id_checked",id_checked),("embedded_id_match_rate",id_rate),("liftover_mapped",lifted_n),("liftover_success_rate",lift_rate),("1kg_reference_rows",ref_rows),("candidate_reference_matches",candidate_n),("duplicate_final_1kg_ids",len(dup_refmap)),("kept_reference_compatible",kept),("retention_vs_additive",retention),("kept_MHCexcluded",kept_mhc)]
    for k,v in sorted(reasons.items()): rows.append(("drop_"+k,v))
    for k,v in sorted(allele_status.items()): rows.append(("allele_"+k,v))
    pd.DataFrame(rows,columns=["metric","value"]).to_csv(out/cfg["outputs"]["harmonization_report"],sep="\t",index=False)
    pd.DataFrame([{"raw_build":"GRCh38","target_build":"GRCh37","strategy":cfg["harmonization"]["strategy"],"embedded_id_match_rate":id_rate,"liftover_success_rate":lift_rate,"reference_compatible_retention":retention,"kept_n":kept}]).to_csv(out/cfg["outputs"]["liftover_audit"],sep="\t",index=False)
    write_text(out/cfg["outputs"]["harmonization_pass"],f"status=PASS\nstage=MDV7_SPARK_HARMONIZATION\nassociation_rows_read=YES\nraw_build=GRCh38\ntarget_build=GRCh37\nadditive_rows={add_rows}\nliftover_success_rate={lift_rate:.8f}\nreference_compatible_snps={kept}\nretention_vs_additive={retention:.8f}\nsample_size_mode=PER_VARIANT_OBS_CT\n")
    import shutil
    shutil.rmtree(tmp,ignore_errors=True)
    print("MDV7_HARMONIZATION_PASS",kept)
if __name__=="__main__": main()
