#!/usr/bin/env bash
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CMD="${1:-status}"
ARG="${2:-}"
PY="${MAR04_PYTHON:-$(command -v python3)}"

test -n "$PY" && test -x "$PY" || { echo "HOLD: python3 not found"; exit 20; }

case "$CMD" in
  preflight|freeze|analysis|postflight|pack-review|status)
    exec "$PY" "$ROOT/workflow/mar04.py" "$CMD"
    ;;
  release)
    [ "$ARG" = "RELEASE" ] || { echo "Usage: bash run_mar04.sh release RELEASE"; exit 2; }
    exec "$PY" "$ROOT/workflow/mar04.py" release
    ;;
  all)
    "$PY" "$ROOT/workflow/mar04.py" preflight || exit $?
    "$PY" "$ROOT/workflow/mar04.py" freeze || exit $?
    OUT="/home/h3021/chapter3/12_molecular_autism_revision/04_effect_boundary/MAR04_effect_boundary_v1"
    if [ ! -f "$OUT/MAR04_STOPA_RELEASE.txt" ]; then
      echo "HOLD_FOR_MANUAL_GATE"
      echo "Review MAR04_preanalysis_lock.json, then run:"
      echo "  bash run_mar04.sh release RELEASE"
      echo "After release, rerun:"
      echo "  bash run_mar04.sh all-post"
      exit 0
    fi
    "$PY" "$ROOT/workflow/mar04.py" analysis || exit $?
    "$PY" "$ROOT/workflow/mar04.py" postflight || exit $?
    "$PY" "$ROOT/workflow/mar04.py" pack-review
    ;;
  all-post)
    "$PY" "$ROOT/workflow/mar04.py" analysis || exit $?
    "$PY" "$ROOT/workflow/mar04.py" postflight || exit $?
    "$PY" "$ROOT/workflow/mar04.py" pack-review
    ;;
  *)
    echo "Usage: bash run_mar04.sh {preflight|freeze|release RELEASE|analysis|postflight|pack-review|status|all|all-post}"
    exit 2
    ;;
esac
