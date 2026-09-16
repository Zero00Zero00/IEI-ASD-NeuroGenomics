#!/usr/bin/env python3
from common import *
import tarfile
if not pass_file(OUT/'MAR02_PASS.txt'): raise SystemExit('HOLD: MAR02 not PASS')
H=MA/'90_handoff'; H.mkdir(parents=True,exist_ok=True); stamp=dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ'); tar=H/f'MAR02_REVIEW_{stamp}.tar.gz'
keep=['MAR02_raw26_redundancy_map.tsv','MAR02_raw26_coclustering.tsv.gz','MAR02_raw26_units.tsv','MAR02_raw26_domain_membership.tsv','MAR02_raw26_standalones.tsv','MAR02_raw26_vs_legacy_mdv5.tsv','MAR02_raw26_vs_legacy_mdv5_summary.tsv','MAR02_raw26_threshold_sensitivity.tsv','MAR02_COMMON_HYPOTHESIS_LOCK.json','MAR02_preflight_report.tsv','MAR02_preanalysis_lock.json','MAR02_postflight_checks.tsv','MAR02_validation_summary.txt','MAR02_checksums.sha256','MAR02_PASS.txt']
with tarfile.open(tar,'w:gz') as tf:
    for n in keep:
        p=OUT/n
        if p.is_file(): tf.add(p,arcname=n)
sha=Path(str(tar)+'.sha256'); sha.write_text(f'{sha256_file(tar)}  {tar.name}\n')
print(f'REVIEW_PACKAGE={tar}\nSHA256={sha}\nPACK_REVIEW=PASS')
