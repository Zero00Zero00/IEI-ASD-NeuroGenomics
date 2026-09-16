#!/usr/bin/env python3
from common import *
require_release(); cfg=WORK/'MAR02_adapter_config.yaml'; LOGS.mkdir(parents=True,exist_ok=True)
for script,name in [('mdv5_02_jaccard.py','01_jaccard'),('mdv5_03_redundancy.py','02_redundancy'),('mdv5_04_leiden.py','03_leiden'),('mdv5_05_domain_freeze.py','04_domain_freeze')]:
    print(f'RUN {name} ...',flush=True); run_vendor(script,cfg,LOGS/f'{name}.log'); print(f'PASS {name}',flush=True)
atomic_text(OUT/'MAR02_ANALYSIS_COMPLETE.txt',f'status=PASS\nstage=MAR02_ANALYSIS\ntimestamp_utc={utc()}\nassociation_results_read=NO\ngwas_inputs_used=NO')
print('MAR02_ANALYSIS=PASS')
