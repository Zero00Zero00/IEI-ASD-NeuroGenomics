#!/usr/bin/env python3
from common import *
if not pass_file(OUT/'MAR02_CANONICALIZE_PASS.txt'): raise SystemExit('HOLD: canonical outputs missing')

def pathway_set(p, key='pathway_key'):
    d=read_tsv(p); return set(d[key].astype(str))
def sethash(s): return hash_gene_set(sorted(s))
def unit_sets(p, flag, unitcol):
    d=read_tsv(p); gcol='HGNC_id' if 'HGNC_id' in d.columns else 'gene'; out={}
    for u,x in d.groupby(unitcol):
        mask=pd.to_numeric(x[flag],errors='coerce').fillna(0).astype(int)==1
        out[str(u)]=set(x.loc[mask,gcol].astype(str))
    return out
def best_match(curset, legacy):
    if not legacy:return ('',set(),0.0)
    scored=[]
    for uid,s in legacy.items():
        j=len(curset&s)/len(curset|s) if curset|s else 1.0; scored.append((j,uid,s))
    j,uid,s=max(scored,key=lambda x:(x[0],x[1])); return uid,s,j
rows=[]; summary=[]
# Pathway-level sets.
cur_t=set(canon_pathways(read_tsv(WP01/'WP01_raw26_tier1.tsv'))['pathway_key']); leg_t=pathway_set(LEGACY/'MDV5_tier1_annotated.tsv')
cur_r=pathway_set(OUT/'MAR02_raw26_representatives.tsv'); leg_r=pathway_set(LEGACY/'MDV5_representative_pathways.tsv')
cur_s=pathway_set(OUT/'MAR02_raw26_standalones.tsv'); leg_s=pathway_set(LEGACY/'MDV5_standalone_pathways.tsv')
for typ,a,b in [('TIER1_PATHWAY_SET',cur_t,leg_t),('REPRESENTATIVE_PATHWAY_SET',cur_r,leg_r),('STANDALONE_PATHWAY_SET',cur_s,leg_s)]:
    rows.append([typ,'COLLECTION','COLLECTION',sethash(a),sethash(b),len(a),len(b),int(a==b),len(a-b),';'.join(sorted(a-b)),len(b-a),';'.join(sorted(b-a)),len(a&b)/len(a|b) if a|b else 1.0])
# Higher-order unit collection by Expanded / Compact exact memberships.
curdm=OUT/'MAR02_raw26_domain_membership.tsv'; legdm=LEGACY/'06_domain_membership.tsv'
for flag,typ in [('Expanded_flag','HIGHER_ORDER_EXPANDED'),('Compact_flag','HIGHER_ORDER_COMPACT')]:
    cur=unit_sets(curdm,flag,'unit_id'); leg=unit_sets(legdm,flag,'domain_id')
    matched=[]
    for uid,s in sorted(cur.items()):
        exact=[lid for lid,ls in leg.items() if ls==s]
        if exact:
            lid=sorted(exact)[0]; ls=leg[lid]; j=1.0
        else:
            lid,ls,j=best_match(s,leg)
        matched.append((uid,lid,s,ls,j))
        rows.append([typ,uid,lid,sethash(s),sethash(ls),len(s),len(ls),int(s==ls),len(s-ls),';'.join(sorted(s-ls)),len(ls-s),';'.join(sorted(ls-s)),j])
    coll_exact=(len(cur)==len(leg) and sorted(sethash(s) for s in cur.values())==sorted(sethash(s) for s in leg.values()))
    summary.append([typ,len(cur),len(leg),int(coll_exact)])
cols=['entity','raw26_id','legacy_id','raw26_hash','legacy_hash','raw26_n','legacy_n','exact_same','added_n','added_members','removed_n','removed_members','jaccard']
write_tsv(pd.DataFrame(rows,columns=cols),OUT/'MAR02_raw26_vs_legacy_mdv5.tsv')
summary += [['TIER1_PATHWAY_SET',len(cur_t),len(leg_t),int(cur_t==leg_t)],['REPRESENTATIVE_PATHWAY_SET',len(cur_r),len(leg_r),int(cur_r==leg_r)],['STANDALONE_PATHWAY_SET',len(cur_s),len(leg_s),int(cur_s==leg_s)]]
sdf=pd.DataFrame(summary,columns=['component','raw26_n','legacy_n','exact_same']); allsame=bool((sdf.exact_same==1).all()); sdf.loc[len(sdf)]=['ALL_DOWNSTREAM_MEMBERSHIP',None,None,int(allsame)]; write_tsv(sdf,OUT/'MAR02_raw26_vs_legacy_mdv5_summary.tsv')
atomic_text(OUT/'MAR02_LEGACY_COMPARE_PASS.txt',f'status=PASS\nall_downstream_membership_exact_same={str(allsame).upper()}')
print(f'LEGACY_COMPARE=PASS\nALL_DOWNSTREAM_EXACT_SAME={str(allsame).upper()}')
