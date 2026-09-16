#!/usr/bin/env bash
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MA="${MA_ROOT:-/home/h3021/chapter3/12_molecular_autism_revision}"
PY="${MAR08_PYTHON:-$MA/envs/mar02_mdv5_legacy/bin/python}"
CMD="${1:-status}"; ARG="${2:-}"
test -x "$PY" || { echo "HOLD: validated Python missing: $PY"; exit 20; }
run(){ "$PY" "$@"; }
case "$CMD" in
 preflight) run "$ROOT/workflow/preflight.py" ;;
 freeze) run "$ROOT/workflow/freeze.py" ;;
 release) [ "$ARG" = "RELEASE" ] || { echo "Usage: bash run_mar08_r1p1.sh release RELEASE"; exit 2; }; run "$ROOT/workflow/release.py" RELEASE ;;
 analysis) run "$ROOT/workflow/analysis_abc.py" && run "$ROOT/workflow/analysis_influence.py" ;;
 sensitivity) run "$ROOT/workflow/recurrence.py" ;;
 postflight) run "$ROOT/workflow/postflight.py" ;;
 pack-review) run "$ROOT/workflow/pack_review.py" ;;
 status) run "$ROOT/workflow/status.py" ;;
 all)
   run "$ROOT/workflow/preflight.py" || exit $?
   run "$ROOT/workflow/freeze.py" || exit $?
   OUT="$MA/08_driver_influence/MAR08_failure_localization_influence_v1p1"
   if [ ! -f "$OUT/MA_MAR08_STOPA_RELEASE.txt" ]; then
     echo "HOLD_FOR_MANUAL_GATE"
     echo "Review: $OUT/MAR08_preanalysis_lock.json"
     echo "Review: $OUT/MAR08_FROZEN_SCENARIO_MANIFEST.tsv"
     echo "Then: bash run_mar08_r1p1.sh release RELEASE"
     exit 0
   fi ;;
 all-post)
   run "$ROOT/workflow/analysis_abc.py" || exit $?
   run "$ROOT/workflow/analysis_influence.py" || exit $?
   run "$ROOT/workflow/recurrence.py" || exit $?
   run "$ROOT/workflow/postflight.py" || exit $?
   run "$ROOT/workflow/pack_review.py" ;;
 *) echo "Usage: bash run_mar08_r1p1.sh {preflight|freeze|release RELEASE|analysis|sensitivity|postflight|pack-review|status|all|all-post}"; exit 2 ;;
esac
