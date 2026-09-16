
from common import *
import sys

lock=OUT/"MAR03_preanalysis_lock.json"
if not lock.is_file(): raise SystemExit("HOLD: freeze first")
if len(sys.argv)<2 or sys.argv[1]!="RELEASE":
    raise SystemExit("Usage: release.py RELEASE")

token={
 "decision":"PASS",
 "track":"MOLECULAR_AUTISM_REVISION",
 "stage":"MAR03",
 "gate":"A",
 "reviewer":"MANUAL_EXPERT_REVIEW",
 "timestamp_utc":utc(),
 "preanalysis_lock_sha256":sha256_file(lock),
 "notes":"7 raw26 higher-order Expanded units + 13 standalones frozen; real MDV6/MDV7 GSA and MDV8 matching interfaces locked before association score parsing."
}
(OUT/"MAR03_STOPA_RELEASE.txt").write_text("\n".join(f"{k}={v}" for k,v in token.items())+"\n",encoding="utf-8")
print("MAR03_GATE_A=PASS")
