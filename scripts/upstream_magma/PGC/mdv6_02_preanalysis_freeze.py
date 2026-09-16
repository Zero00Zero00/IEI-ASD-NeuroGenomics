#!/usr/bin/env python3
import argparse, json, os
from pathlib import Path
import pandas as pd
from mdv6_common import *

def setline(name, ids): return name+"\t"+"\t".join(sorted(map(str,ids)))+"\n"

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); a=ap.parse_args()
    cfg=load_cfg(a.config); root=Path(cfg["project_root"]); out=out_dir(cfg); prov=prov_dir(cfg)
    require_file(out/cfg["outputs"]["preflight_pass"])
    u=pd.read_csv(root_path(cfg,cfg["upstream"]["gene_universe"]),sep="\t")
    d=pd.read_csv(root_path(cfg,cfg["upstream"]["mdv5_domain_membership"]),sep="\t")
    lab=pd.read_csv(root_path(cfg,cfg["upstream"]["mdv5_biological_labels"]),sep="\t")
    # MAGMA IDs: stable numeric HGNC suffix, with explicit map.
    u=u.copy(); u["MAGMA_GENE"]=u["HGNC_id"].map(hgnc_numeric)
    if u["MAGMA_GENE"].duplicated().any(): raise RuntimeError("MAGMA_GENE not unique")
    idmap=u[["MAGMA_GENE","HGNC_id","HGNC_symbol","gencode_gene_id","chr","start","end","strand","CoreSeed_flag"]].copy()
    idmap.to_csv(out/cfg["outputs"]["gene_id_map"],sep="\t",index=False)
    # Gene locations: chromosomes 1-23 (X=23); Y/MT are audit-excluded because standard MAGMA EUR ref is not guaranteed to support them.
    loc=[]; excluded=[]
    for _,r in u.iterrows():
        ch=norm_chr(r["chr"])
        if ch is None or ch<1 or ch>23:
            excluded.append((r["HGNC_id"],r["HGNC_symbol"],r["chr"])); continue
        loc.append((r["MAGMA_GENE"],ch,int(r["start"]),int(r["end"]),r["strand"],r["HGNC_symbol"]))
    with open(out/cfg["outputs"]["gene_locations"],"w",encoding="utf-8") as f:
        for x in loc: f.write("\t".join(map(str,x))+"\n")
    # gene sets
    h2m=dict(zip(u.HGNC_id,u.MAGMA_GENE)); sym2m=dict(zip(u.HGNC_symbol,u.MAGMA_GENE))
    units=cfg["expected"]["unit_ids"]
    sets={}; compact={}; alltier={}; removed={}
    core=set(u.loc[pd.to_numeric(u.CoreSeed_flag,errors="coerce").fillna(0).eq(1),"HGNC_id"])
    for uid in units:
        sub=d[d.domain_id==uid]
        sets[uid]=set(sub.loc[sub.Expanded_flag==1,"HGNC_id"])
        compact[uid]=set(sub.loc[sub.Compact_flag==1,"HGNC_id"])
        alltier[uid]=set(sub.loc[sub.AllTier1_flag==1,"HGNC_id"])
        removed[uid]=sets[uid]-core
        if len(sets[uid])!=cfg["expected"]["expanded_counts"][uid]: raise RuntimeError(f"{uid} Expanded count drift")
        if len(compact[uid])!=cfg["expected"]["compact_counts"][uid]: raise RuntimeError(f"{uid} Compact count drift")
    # overlap audit and hard nested check
    rows=[]
    for i,a1 in enumerate(units):
        for b1 in units[i+1:]:
            for typ,dd in [("Expanded",sets),("Compact",compact)]:
                A=dd[a1]; B=dd[b1]; inter=A&B; union=A|B
                rows.append({"set_type":typ,"unit_a":a1,"unit_b":b1,"n_a":len(A),"n_b":len(B),"intersection_n":len(inter),"jaccard":len(inter)/len(union) if union else 0,"a_subset_b":A<=B,"b_subset_a":B<=A})
    ov=pd.DataFrame(rows); ov.to_csv(out/cfg["outputs"]["unit_overlap_audit"],sep="\t",index=False)
    if not sets["M03"] <= sets["M04"]: raise RuntimeError("Frozen expected nesting M03 Expanded subset M04 failed")
    if not compact["M03"] <= compact["M04"]: raise RuntimeError("Frozen expected nesting M03 Compact subset M04 failed")
    nested={uid:set(v) for uid,v in sets.items() if uid!="M04"}
    nested["M04_EXCLUSIVE"]=sets["M04"]-sets["M03"]
    # output set files using numeric HGNC MAGMA IDs
    def write_sets(path, dd):
        with open(path,"w",encoding="utf-8") as f:
            for name in dd:
                ids=[h2m[x] for x in dd[name] if x in h2m]
                f.write(setline(name,ids))
    write_sets(out/cfg["outputs"]["primary_sets"],sets)
    write_sets(out/cfg["outputs"]["compact_sets"],compact)
    write_sets(out/cfg["outputs"]["remove_coreseed_sets"],removed)
    write_sets(out/cfg["outputs"]["all_tier1_sets"],alltier)
    write_sets(out/cfg["outputs"]["nested_sets"],nested)
    with open(out/cfg["outputs"]["coreseed_set"],"w",encoding="utf-8") as f: f.write(setline("CoreSeed",[h2m[x] for x in core]))
    # combined primary set file includes CoreSeed first then 7 units, created for MAGMA marginal run.
    combined=out/"MDV6_primary_with_CoreSeed.sets"
    with open(combined,"w",encoding="utf-8") as dst:
        dst.write((out/cfg["outputs"]["coreseed_set"]).read_text())
        dst.write((out/cfg["outputs"]["primary_sets"]).read_text())
    # summary with biological labels
    lmap=lab.set_index("unit_id").to_dict("index")
    srows=[]
    for uid in units:
        srows.append({"unit_id":uid,"biological_label":lmap[uid]["biological_label"],"driver_architecture":lmap[uid]["driver_architecture"],"Expanded_n":len(sets[uid]),"Compact_n":len(compact[uid]),"Expanded_remove_CoreSeed_n":len(removed[uid]),"AllTier1_n":len(alltier[uid]),"nested_note":"M03 is subset of M04" if uid in ["M03","M04"] else ""})
    pd.DataFrame(srows).to_csv(out/cfg["outputs"]["gene_set_summary"],sep="\t",index=False)
    # freeze hashes incl code/config + set definitions; NO PGC association rows read here.
    scripts=sorted((Path(a.config).parent.parent/"workflow"/"scripts").glob("mdv6_*.py"))
    package_root=Path(a.config).parent.parent
    hashes={str(p.relative_to(package_root)):sha256_file(p) for p in [Path(a.config), package_root/"Snakefile.mdv6", *scripts] if p.exists()}
    files=[cfg["outputs"][x] for x in ["gene_id_map","gene_locations","unit_overlap_audit","gene_set_summary","primary_sets","compact_sets","remove_coreseed_sets","all_tier1_sets","nested_sets","coreseed_set"]]+["MDV6_primary_with_CoreSeed.sets"]
    freeze={"version":cfg["version"],"association_results_read":"NO","unit_ids":units,"M03_Expanded_subset_M04":True,"M03_Compact_subset_M04":True,"primary_window_kb":cfg["primary_inference"]["gene_window_kb"],"primary_gene_sets":"Domain-Expanded","Gate_B":cfg["primary_inference"]["gate_b_primary"],"interpretation_lock":cfg["interpretation_lock"],"scientific_files":{fn:sha256_file(out/fn) for fn in files},"code_and_config_sha256":hashes,"excluded_gene_location_records":excluded}
    write_json(out/cfg["outputs"]["preanalysis_lock"],freeze)
    write_text(out/cfg["outputs"]["preanalysis_pass"],"status=PASS\ntechnical_status=PASS\nstage=MDV6_PREANALYSIS_FREEZE\nassociation_results_read=NO\nfrozen_units=7\nprimary_window_kb=10\nprimary_set=Domain-Expanded\nM03_nested_in_M04=YES\n")
    print("MDV6_PREANALYSIS_FREEZE_PASS")
if __name__=="__main__": main()
