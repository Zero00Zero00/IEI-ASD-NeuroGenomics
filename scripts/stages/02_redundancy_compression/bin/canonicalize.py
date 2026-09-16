#!/usr/bin/env python3
from common import *
if not pass_file(OUT/'MAR02_ANALYSIS_COMPLETE.txt'): raise SystemExit('HOLD: analysis incomplete')

def rn(df):
    mp={'CoreSeed_driver_n':'raw26_driver_n','CoreSeed_drivers':'raw26_drivers','CoreSeed_flag':'raw26_anchor_flag','DomainCore_flag':'raw26_domain_core_flag','DomainCore_gene_n':'raw26_domain_core_gene_n','unique_CoreSeed_driver_n':'unique_raw26_driver_n'}
    return df.rename(columns={k:v for k,v in mp.items() if k in df.columns})

rmap=rn(read_tsv(WORK/'internal_redundancy_map.tsv')); write_tsv(rmap,OUT/'MAR02_raw26_redundancy_map.tsv')
reps=rn(read_tsv(WORK/'internal_representatives.tsv')); write_tsv(reps,OUT/'MAR02_raw26_representatives.tsv')
cp=read_tsv(WORK/'internal_coclustering_pairs.tsv.gz'); write_tsv(cp,OUT/'MAR02_raw26_coclustering.tsv.gz')
ptu=read_tsv(WORK/'internal_pathway_to_unit.tsv'); write_tsv(ptu,OUT/'MAR02_raw26_pathway_to_unit.tsv')
dm=rn(read_tsv(WORK/'internal_domain_membership.tsv')); dm=dm.rename(columns={'domain_id':'unit_id','domain_label':'biological_label_placeholder'}); write_tsv(dm,OUT/'MAR02_raw26_domain_membership.tsv')
ds=rn(read_tsv(WORK/'internal_domain_summary.tsv')); ds=ds.rename(columns={'domain_id':'unit_id','domain_label':'biological_label_placeholder','representative_pathway_n':'representative_count','original_tier1_pathway_n':'original_count','Expanded_gene_n':'Expanded_n','Compact_gene_n':'Compact_n','AllTier1_gene_n':'AllTier1_n'}); write_tsv(ds,OUT/'MAR02_raw26_units.tsv')
st=rn(read_tsv(WORK/'internal_standalones.tsv')).sort_values('pathway_key',kind='mergesort').reset_index(drop=True); st.insert(0,'standalone_id',[f'SA{i:03d}' for i in range(1,len(st)+1)]); write_tsv(st,OUT/'MAR02_raw26_standalones.tsv')
ld=rn(read_tsv(WORK/'internal_locus_dominance.tsv')).rename(columns={'domain_id':'unit_id'}); write_tsv(ld,OUT/'MAR02_raw26_locus_dominance.tsv')
for src,dst in [('internal_threshold_sensitivity.tsv','MAR02_raw26_threshold_sensitivity.tsv'),('internal_network_edges.tsv','MAR02_raw26_network_edges.tsv'),('internal_network_edges_sensitivity.tsv.gz','MAR02_raw26_network_edges_sensitivity.tsv.gz'),('internal_jaccard_pairs.tsv.gz','MAR02_raw26_jaccard_pairs.tsv.gz')]:
    write_tsv(read_tsv(WORK/src),OUT/dst)
atomic_text(OUT/'MAR02_CANONICALIZE_PASS.txt',f'status=PASS\ntimestamp_utc={utc()}')
print(f'CANONICALIZE=PASS\nUNITS={len(ds)}\nSTANDALONES={len(st)}\nREPRESENTATIVES={len(reps)}')
