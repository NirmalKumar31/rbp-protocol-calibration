#!/usr/bin/env bash
# Rebuild COLUMNS.csv and PROVENANCE.csv, in the only order that converges.
#
# The two generators are circular. provenance.py records a sha256 of COLUMNS.csv, and
# column_dictionary.py reads producing_script out of PROVENANCE.csv. So a newly committed table
# comes out of the first pass with no producer, and provenance then hashes that incomplete file.
# Both --check modes still pass locally, because the two files agree with EACH OTHER; only a
# fresh build disagrees with both, which means CI finds it and you do not.
#
# Two passes fix it: the producer map is settled by the end of the first provenance run.
# Run this instead of calling either script by hand, and run it LAST, after every other
# generator, because anything regenerated afterwards leaves the manifest stale again.
set -euo pipefail
cd "$(dirname "$0")/.."
PY="${PY:-python3}"
export PYTHONPATH="$PWD/src"

for pass in 1 2; do
  "$PY" scripts/column_dictionary.py >/dev/null
  "$PY" scripts/provenance.py >/dev/null
done

"$PY" scripts/column_dictionary.py --check
"$PY" scripts/provenance.py --check
