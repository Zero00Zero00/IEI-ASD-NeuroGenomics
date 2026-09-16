#!/usr/bin/env bash
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MA="${MA_ROOT:-/home/h3021/chapter3/12_molecular_autism_revision}"
PY="${MAR06_PYTHON:-$MA/envs/mar02_mdv5_legacy/bin/python}"
CMD="${1:-status}"
ARG="${2:-}"
test -x "$PY" || { echo "HOLD: validated Python missing: $PY"; exit 20; }
run(){ "$PY" "$@"; }

case "$CMD" in
  discover) run "$ROOT/workflow/discover_raw26.py" ;;
  preflight) run "$ROOT/workflow/preflight.py" ;;
  freeze) run "$ROOT/workflow/freeze.py" ;;
  release)
    [ "$ARG" = "RELEASE" ] || { echo "Usage: bash run_mar06.sh release RELEASE"; exit 2; }
    run "$ROOT/workflow/release.py" RELEASE ;;
  raw26) run "$ROOT/workflow/audit_raw26.py" ;;
  common) run "$ROOT/workflow/audit_common.py" ;;
  postflight) run "$ROOT/workflow/postflight.py" ;;
  pack-review) run "$ROOT/workflow/pack_review.py" ;;
  status) run "$ROOT/workflow/status.py" ;;
  all)
    run "$ROOT/workflow/discover_raw26.py" || exit $?
    run "$ROOT/workflow/preflight.py" || exit $?
    run "$ROOT/workflow/freeze.py" || exit $?
    OUT="$MA/06_matched_null_audit/MAR06_matched_null_audit_v1.0R2"
    if [ ! -f "$OUT/MAR06_STOPA_RELEASE.txt" ]; then
      echo "HOLD_FOR_MANUAL_GATE"
      echo "Review MAR06_preanalysis_lock.json, then:"
      echo "  bash run_mar06.sh release RELEASE"
      echo "After release:"
      echo "  bash run_mar06.sh all-post"
      exit 0
    fi
    ;;
  all-post)
    run "$ROOT/workflow/audit_raw26.py" || exit $?
    run "$ROOT/workflow/audit_common.py" || exit $?
    run "$ROOT/workflow/postflight.py" || exit $?
    run "$ROOT/workflow/pack_review.py"
    ;;
  *)
    echo "Usage: bash run_mar06.sh {discover|preflight|freeze|release RELEASE|raw26|common|postflight|pack-review|status|all|all-post}"
    exit 2 ;;
esac
