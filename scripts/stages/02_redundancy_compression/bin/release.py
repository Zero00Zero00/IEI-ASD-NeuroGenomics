#!/usr/bin/env python3
from common import *
arg=sys.argv[1] if len(sys.argv)>1 else ''
lock=OUT/'MAR02_preanalysis_lock.json'
if not lock.is_file(): raise SystemExit('HOLD: run freeze first')
if arg!='RELEASE':
    print('MANUAL GATE A REVIEW REQUIRED')
    print('Confirm: raw26/WP01 inputs only; no PGC/SPARK gene-level results; fixed rules; output namespace MAR02_raw26_compression_v1.')
    print('If accepted, run: bash run_mar02.sh release RELEASE')
    raise SystemExit(4)
atomic_text(OUT/'MAR02_STOPA_RELEASE.txt',f'decision=PASS\ntrack=MOLECULAR_AUTISM_REVISION\nstage=MAR02\ngate=A\nreviewer=MANUAL_EXPERT_REVIEW\ntimestamp_utc={utc()}\npreanalysis_lock_sha256={sha256_file(lock)}\nnotes=Reviewed real-path WP01 raw26 inputs, legacy MDV5 parameter authority, GWAS-blind firewall, and fixed MAR02 rules.')
print('MAR02_GATE_A=PASS')
