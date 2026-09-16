#!/usr/bin/env bash
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MA="${MA_ROOT:-/home/h3021/chapter3/12_molecular_autism_revision}"
PY="${MAR07_PYTHON:-$MA/envs/mar02_mdv5_legacy/bin/python}"
CMD="${1:-status}"; ARG="${2:-}"
test -x "$PY" || { echo "HOLD: validated Python missing: $PY"; exit 20; }
run(){ "$PY" "$@"; }

case "$CMD" in
  preflight) run "$ROOT/workflow/preflight.py" ;;
  freeze) run "$ROOT/workflow/freeze.py" ;;
  release)
    [ "$ARG" = "RELEASE" ] || { echo "Usage: bash run_mar07_r4p1.sh release RELEASE"; exit 2; }
    run "$ROOT/workflow/release.py" RELEASE ;;
  observed) run "$ROOT/workflow/import_observed.py" ;;
  certify) run "$ROOT/workflow/certify.py" ;;
  pilot) run "$ROOT/workflow/calibration.py" pilot ;;
  calibrate) run "$ROOT/workflow/calibration.py" calibrate ;;
  postflight) run "$ROOT/workflow/postflight.py" ;;
  pack-review) run "$ROOT/workflow/pack_review.py" ;;
  status) run "$ROOT/workflow/status.py" ;;
  all)
    run "$ROOT/workflow/preflight.py" || exit $?
    run "$ROOT/workflow/freeze.py" || exit $?
    OUT="$MA/07_chain_calibration/MAR07_chain_calibration_v1.0R4p1"
    if [ ! -f "$OUT/MAR07_STOPA_RELEASE.txt" ]; then
      echo "HOLD_FOR_MANUAL_GATE"
      echo "Review MAR07_preanalysis_lock.json, then:"
      echo "  bash run_mar07_r4p1.sh release RELEASE"
      exit 0
    fi ;;
  all-post)
    run "$ROOT/workflow/import_observed.py" || exit $?
    run "$ROOT/workflow/certify.py" || exit $?
    run "$ROOT/workflow/calibration.py" pilot || exit $?
    run "$ROOT/workflow/calibration.py" calibrate || exit $?
    run "$ROOT/workflow/postflight.py" || exit $?
    run "$ROOT/workflow/pack_review.py" ;;
  *)
    echo "Usage: bash run_mar07_r4p1.sh {preflight|freeze|release RELEASE|observed|certify|pilot|calibrate|postflight|pack-review|status|all|all-post}"
    exit 2 ;;
esac
