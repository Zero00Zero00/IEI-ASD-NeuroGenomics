#!/usr/bin/env bash
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MA="${MA_ROOT:-/home/h3021/chapter3/12_molecular_autism_revision}"
PY="${MAR05_PYTHON:-$MA/envs/mar02_mdv5_legacy/bin/python}"
CMD="${1:-status}"
ARG="${2:-}"

test -x "$PY" || { echo "HOLD: validated Python not found: $PY"; exit 20; }

case "$CMD" in
  preflight) exec "$PY" "$ROOT/workflow/preflight.py" ;;
  freeze) exec "$PY" "$ROOT/workflow/freeze.py" ;;
  release)
    [ "$ARG" = "RELEASE" ] || { echo "Usage: bash run_mar05.sh release RELEASE"; exit 2; }
    exec "$PY" "$ROOT/workflow/release.py" RELEASE ;;
  analysis) exec "$PY" "$ROOT/workflow/analysis.py" ;;
  postflight) exec "$PY" "$ROOT/workflow/postflight.py" ;;
  pack-review) exec "$PY" "$ROOT/workflow/pack_review.py" ;;
  status) exec "$PY" "$ROOT/workflow/status.py" ;;
  all)
    "$PY" "$ROOT/workflow/preflight.py" || exit $?
    "$PY" "$ROOT/workflow/freeze.py" || exit $?
    OUT="$MA/05_magma_collinearity/MAR05_collinearity_v1"
    if [ ! -f "$OUT/MAR05_STOPA_RELEASE.txt" ]; then
      echo "HOLD_FOR_MANUAL_GATE"
      echo "Review MAR05_preanalysis_lock.json, then:"
      echo "  bash run_mar05.sh release RELEASE"
      echo "After release:"
      echo "  bash run_mar05.sh all-post"
      exit 0
    fi
    "$PY" "$ROOT/workflow/analysis.py" || exit $?
    "$PY" "$ROOT/workflow/postflight.py" || exit $?
    "$PY" "$ROOT/workflow/pack_review.py"
    ;;
  all-post)
    "$PY" "$ROOT/workflow/analysis.py" || exit $?
    "$PY" "$ROOT/workflow/postflight.py" || exit $?
    "$PY" "$ROOT/workflow/pack_review.py"
    ;;
  *)
    echo "Usage: bash run_mar05.sh {preflight|freeze|release RELEASE|analysis|postflight|pack-review|status|all|all-post}"
    exit 2 ;;
esac
