#!/usr/bin/env bash
# Run the test suite against the file set the CONTAINER actually has, before paying for a build.
#
#   bash scripts/check_image_tree.sh
#
# WHY THIS EXISTS. Two tests of ours assumed a developer checkout and failed every Cloud Build
# of the GPU image for weeks: one shelled out to `git ls-files` at collection time, and one
# asserted that manuscript figures exist. Neither is visible locally, because locally git and
# the manuscript are both present. The consequence was not a red build anyone watched -- it was
# a STALE PUBLISHED IMAGE, so a Batch job months later died on an arm the image had never heard
# of, and the error pointed at the arm rather than at the build.
#
# The Dockerfile copies config/, src/, scripts/, tests/, docker/bake_weights.py and
# pyproject.toml, and nothing else: no .git, no .gitignore, no cloud/, no manuscript/, no
# results/. This mirrors exactly that set into a temporary tree, hides git from PATH, and runs
# the suite there. It costs under a minute and it is the difference between finding this in a
# shell and finding it in a build log an hour later.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1

PY="${PY:-python3}"
SIM=$(mktemp -d)
trap 'rm -rf "$SIM"' EXIT

mkdir -p "$SIM/docker"
cp -R config src scripts tests "$SIM/"
cp docker/bake_weights.py "$SIM/docker/"
cp pyproject.toml "$SIM/"

echo "running the suite against the image's file set (no git, no .gitignore, no cloud/,"
echo "no manuscript/, no results/) in $SIM"
( cd "$SIM" && env PATH="/nonexistent" PYTHONPATH="$SIM/src" \
    "$PY" -m pytest tests -q -p no:cacheprovider ) || {
  echo
  echo "FAIL: the suite does not pass against the container's file set, so the image build"
  echo "will fail and the PUBLISHED IMAGE WILL SILENTLY STAY STALE. Fix the test to skip"
  echo "or degrade where the file it needs is absent, rather than removing it."
  exit 1
}
echo "OK: the image build's test step will pass."
