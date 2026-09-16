#!/usr/bin/env bash
set -euo pipefail
MA="${MA_ROOT:-/home/h3021/chapter3/12_molecular_autism_revision}"
PY="${MAR07_PYTHON:-$MA/envs/mar02_mdv5_legacy/bin/python}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec "$PY" "$ROOT/executor/fullchain_executor.py" "$@"
