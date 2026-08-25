#!/usr/bin/env bash
# The gate: is a dataset built in the cloud byte-identical to the one built on the laptop?
#
#   bash cloud/verify_prep.sh                 # every cloud dataset that has a local twin
#   bash cloud/verify_prep.sh QKI             # just one protein
#
# Compares md5, not size and not row count. A container that preprocessed differently would
# differ in a handful of windows out of millions and every summary statistic would agree.
#
# The md5 is read with `gcloud storage hash --hex` rather than parsed out of `ls -L`. An
# earlier hand-rolled parser silently returned "?" for every object and reported three
# spurious mismatches; the lesson was to suspect the verification before the data.
# Project and buckets come from the environment, never from a literal. A hardcoded id is
# how a pipeline ends up only running on its author's account. Override with:
#   export GOOGLE_CLOUD_PROJECT=your-project
PROJECT_ID="${GOOGLE_CLOUD_PROJECT:-$(.venv/bin/python -c 'import sys;sys.path.insert(0,"src");from rbp.utils import cloud;print(cloud.project())')}"
DERIVED="${DERIVED_BUCKET:-${PROJECT_ID}-derived}"
RAW="${RAW_BUCKET:-${PROJECT_ID}-raw}"

set -uo pipefail

DERIVED=${DERIVED}
ONLY=${1:-}

# A case, not `declare -A`. macOS ships bash 3.2, where associative arrays do not exist and
# the subscript silently resolves to nothing -- which under `set -u` aborts the whole script
# before a single file is compared.
local_dir() {
  case "$1" in
    gc)    echo data/processed ;;
    dinuc) echo data/processed_dinucmatch ;;
    *)     echo "" ;;
  esac
}

pass=0; fail=0; skip=0
printf '%-10s %-7s %-7s %-8s %s\n' protein cell arm result md5
for obj in $(gcloud storage ls "gs://${DERIVED}/processed/**/dataset.tsv" 2>/dev/null); do
  rel=${obj#gs://${DERIVED}/processed/}
  arm=${rel%%/*}; rest=${rel#*/}; cell=${rest%%/*}; rest=${rest#*/}; prot=${rest%%/*}
  [ -n "$ONLY" ] && [ "$prot" != "$ONLY" ] && continue

  lf="$(local_dir "$arm")/${cell}/${prot}/dataset.tsv"
  if [ ! -f "$lf" ]; then
    printf '%-10s %-7s %-7s %-8s %s\n' "$prot" "$cell" "$arm" SKIP "no local twin"
    skip=$((skip + 1)); continue
  fi

  # --hex because the default is base64, and a base64/hex mix-up looks exactly like a real
  # mismatch. Do NOT add --no-user-output-enabled: it suppresses --format output too, so
  # the hash comes back empty and every row reads DIFFER. That is failure #9 a second time,
  # from a different direction -- hence the explicit empty check below rather than trusting
  # the comparison.
  cloud=$(gcloud storage hash "$obj" --hex --format='value(md5_hash)' 2>/dev/null)
  local_md5=$(md5 -q "$lf" 2>/dev/null || md5sum "$lf" | cut -d' ' -f1)

  if [ -z "$cloud" ]; then
    printf '%-10s %-7s %-7s %-8s %s\n' "$prot" "$cell" "$arm" ERROR "no md5 returned"
    fail=$((fail + 1))
  elif [ "$cloud" = "$local_md5" ]; then
    printf '%-10s %-7s %-7s %-8s %s\n' "$prot" "$cell" "$arm" MATCH "$cloud"
    pass=$((pass + 1))
  else
    printf '%-10s %-7s %-7s %-8s %s\n' "$prot" "$cell" "$arm" DIFFER "cloud=$cloud local=$local_md5"
    fail=$((fail + 1))
  fi
done

echo
echo "match $pass   differ/error $fail   no local twin $skip"
[ "$fail" -eq 0 ] && [ "$pass" -gt 0 ]
