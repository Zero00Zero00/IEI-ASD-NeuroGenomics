#!/usr/bin/env python3
from common import *
print(f'MAR02_OUT={OUT}')
for name in ['MAR02_preflight_report.tsv','MAR02_preanalysis_lock.json','MAR02_STOPA_RELEASE.txt','MAR02_ANALYSIS_COMPLETE.txt','MAR02_SENSITIVITY_PASS.txt','MAR02_COMMON_HYPOTHESIS_LOCK.json','MAR02_postflight_checks.tsv','MAR02_checksums.sha256','MAR02_PASS.txt']:
    p=OUT/name; print(('OK   ' if p.is_file() else 'MISS ')+str(p))
if (OUT/'MAR02_PASS.txt').is_file(): print('\n'+(OUT/'MAR02_PASS.txt').read_text())
