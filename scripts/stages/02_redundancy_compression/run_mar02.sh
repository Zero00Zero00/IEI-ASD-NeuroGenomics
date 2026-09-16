#!/usr/bin/env bash
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="$ROOT/bin:${PYTHONPATH:-}"
CMD="${1:-help}"
runpy(){ python3 "$ROOT/bin/$1" "${@:2}"; }
case "$CMD" in
  preflight) runpy preflight.py ;;
  freeze) runpy freeze.py ;;
  release) runpy release.py "${2:-}" ;;
  analysis) runpy analysis.py && runpy canonicalize.py && runpy compare_legacy.py ;;
  sensitivity) runpy sensitivity.py ;;
  postflight) runpy postflight.py ;;
  status) runpy status.py ;;
  pack-review) runpy pack_review.py ;;
  all)
    runpy preflight.py || exit $?
    runpy freeze.py || exit $?
    if [ ! -f /home/h3021/chapter3/12_molecular_autism_revision/02_raw26_domains/MAR02_raw26_compression_v1/MAR02_STOPA_RELEASE.txt ]; then
      echo "HOLD_FOR_MANUAL_GATE: inspect MAR02_preanalysis_lock.json, then run: bash run_mar02.sh release RELEASE"
      exit 4
    fi
    runpy analysis.py || exit $?
    runpy canonicalize.py || exit $?
    runpy compare_legacy.py || exit $?
    runpy sensitivity.py || exit $?
    runpy postflight.py || exit $?
    runpy pack_review.py
    ;;
  help|-h|--help)
    cat <<'EOF'
MAR02 raw26 compression — one-click staged runner

Commands:
  bash run_mar02.sh preflight
  bash run_mar02.sh freeze
  bash run_mar02.sh release RELEASE   # manual Gate A after reviewing preanalysis lock
  bash run_mar02.sh analysis
  bash run_mar02.sh sensitivity
  bash run_mar02.sh postflight
  bash run_mar02.sh status
  bash run_mar02.sh pack-review
  bash run_mar02.sh all               # stops at manual gate on first run; rerun after release

Recommended:
  1) bash run_mar02.sh all
  2) review 02_raw26_domains/MAR02_raw26_compression_v1/MAR02_preanalysis_lock.json
  3) bash run_mar02.sh release RELEASE
  4) bash run_mar02.sh all
EOF
    ;;
  *) echo "Unknown command: $CMD"; exit 2;;
esac
