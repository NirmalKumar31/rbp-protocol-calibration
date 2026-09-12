#!/usr/bin/env bash
# Everything CI would run for the dashboard, in one place.
#
# Wired into .github/workflows/ci.yml as the `dashboard` job. It exists as a script rather than
# as inline YAML so the same checks run locally and in CI from one definition.
#
# It is a separate job because dashboard/requirements.txt installs Streamlit and Plotly, which
# the study environment does not have and must not acquire: a Plotly bump cannot be allowed to
# change a published number. The `test` and `full-suite` jobs run `pytest tests`, which never
# reaches dashboard/tests, so before this job the dashboard suite ran nowhere.
set -euo pipefail
cd "$(dirname "$0")/.."

PY="${PY:-python}"
fail() { echo "FAILED: $1" >&2; exit 1; }

echo "=== ruff"
"$PY" -m ruff check dashboard || fail "ruff"

echo "=== dashboard tests"
"$PY" -m pytest dashboard/tests -q || fail "pytest"

echo "=== every view renders at every filter width"
"$PY" - <<'PYEOF' || fail "view rendering"
import sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "dashboard")
import importlib
app = importlib.import_module("app")
from rbp_dashboard import data

frame = data.table("three_arm_per_dataset.csv")
proteins = sorted(frame["protein"].unique())
cases = [
    ("all", {"cells": None, "proteins": None, "size_range": None, "size_column": None}),
    ("one cell line", {"cells": ["K562"], "proteins": [], "size_range": None,
                       "size_column": None}),
    ("one protein", {"cells": None, "proteins": proteins[:1], "size_range": None,
                     "size_column": None}),
    ("empty", {"cells": ["no such cell line"], "proteins": [], "size_range": None,
               "size_column": None}),
]
for label, flt in cases:
    for name, view in app.VIEWS.items():
        view(flt)
    print(f"  {label}: {len(app.VIEWS)} views ok")
PYEOF

echo
echo "ALL DASHBOARD CHECKS PASSED"
