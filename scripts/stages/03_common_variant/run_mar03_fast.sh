#!/usr/bin/env bash
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MA="${MA_ROOT:-/home/h3021/chapter3/12_molecular_autism_revision}"
PY="${MAR03_PYTHON:-$MA/envs/mar02_mdv5_legacy/bin/python}"
CMD="${1:-status}"
ARG="${2:-}"

test -x "$PY" || { echo "HOLD: validated Python not found: $PY"; exit 20; }

run() { "$PY" "$@"; }

case "$CMD" in
  preflight) run "$ROOT/workflow/preflight.py" ;;
  freeze) run "$ROOT/workflow/freeze.py" ;;
  release) run "$ROOT/workflow/release.py" "$ARG" ;;
  analysis) run "$ROOT/workflow/analysis.py" ;;
  matched) run "$ROOT/workflow/matched.py" ;;
  sensitivity) run "$ROOT/workflow/sensitivity.py" ;;
  postflight) run "$ROOT/workflow/postflight.py" ;;
  pack-review) run "$ROOT/workflow/pack_review.py" ;;
  status) run "$ROOT/workflow/status.py" ;;
  all)
    run "$ROOT/workflow/preflight.py" || exit $?
    run "$ROOT/workflow/freeze.py" || exit $?
    echo "HOLD_FOR_MANUAL_GATE"
    echo "Review MAR03_preanalysis_lock.json, then run:"
    echo "  bash run_mar03_fast.sh release RELEASE"
    exit 0
    ;;
  *)
    echo "Usage: bash run_mar03_fast.sh {preflight|freeze|release RELEASE|analysis|matched|sensitivity|postflight|pack-review|status|all}"
    exit 2
    ;;
esac
