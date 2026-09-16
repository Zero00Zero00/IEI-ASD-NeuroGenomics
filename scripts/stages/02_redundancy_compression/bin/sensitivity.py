#!/usr/bin/env python3
from common import *
if not pass_file(OUT/'MAR02_CANONICALIZE_PASS.txt'): raise SystemExit('HOLD: analysis/canonicalization missing')
d=read_tsv(OUT/'MAR02_raw26_threshold_sensitivity.tsv'); vals=sorted(round(float(x),12) for x in d['edge_threshold'].unique()); expected=[0.15,0.20,0.25,0.30]
if vals!=expected: raise SystemExit(f'HOLD: threshold sensitivity mismatch observed={vals} expected={expected}')
write_tsv(d[['edge_threshold','is_primary','edge_n','community_n','ARI_vs_primary_medoid','stable_domain_n','small_module_n','medoid_seed','medoid_loss']],OUT/'MAR02_sensitivity_summary.tsv')
atomic_text(OUT/'MAR02_SENSITIVITY_PASS.txt',f'status=PASS\nthresholds={vals}\nprimary=0.20')
print('MAR02_SENSITIVITY=PASS')
