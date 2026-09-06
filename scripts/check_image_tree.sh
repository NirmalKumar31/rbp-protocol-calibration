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
#
# TWO DEFECTS AN EXTERNAL AUDIT FOUND IN THIS FILE, both of which made it useless:
#
#   1. PY defaulted to a BARE NAME, python3, and the interpreter then ran under
#      `env PATH="/nonexistent"`, which is the whole point of the script. env looks the command
#      up in the PATH it is about to set, so it never found it: `env: python3: No such file or
#      directory`, every time, for everybody. run.sh gates stage 2 on this script, so the
#      documented end-to-end rebuild died before the first paid step. Resolve to an ABSOLUTE
#      path first, while PATH still exists.
#
#   2. It ran one selection -- the whole tests/ tree -- against a simulated CPU image. The CPU
#      Cloud Build deliberately ships without torch and ignores test_models.py and
#      test_train_folds.py; the GPU build runs everything. So the simulation matched NEITHER
#      build: it over-tested the CPU image (collecting torch tests the image cannot import) and
#      under-described the GPU one. Both selections are run below, read from the cloudbuild
#      files rather than copied here, because a third copy of a list is a third thing to drift.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1

PY="$(command -v "${PY:-python3}")"
[ -x "$PY" ] || { echo "FAIL: no usable interpreter (PY=${PY:-python3})"; exit 1; }

SIM=$(mktemp -d)
trap 'rm -rf "$SIM"' EXIT

mkdir -p "$SIM/docker"
cp -R config src scripts tests "$SIM/"
cp docker/bake_weights.py "$SIM/docker/"
cp pyproject.toml "$SIM/"

# The CPU build's ignore list, taken from the build file so the two cannot disagree.
# `mapfile` is bash 4; macOS ships bash 3.2 and this script runs on both, so read a loop.
CPU_IGNORE=()
while IFS= read -r line; do
  [ -n "$line" ] && CPU_IGNORE+=("$line")
done < <(grep -oE -- '--ignore=tests/[A-Za-z0-9_/.]+' docker/cloudbuild.cpu.yaml | sort -u)
[ ${#CPU_IGNORE[@]} -gt 0 ] || { echo "FAIL: no --ignore flags found in cloudbuild.cpu.yaml;"
  echo "either the CPU build stopped excluding the torch tests or this parser broke."; exit 1; }

echo "the image's file set: no git, no .gitignore, no cloud/, no manuscript/, no results/"
echo "simulated in $SIM"

fail=0
run_one() {
  local label=$1; shift
  echo
  echo "--- $label ---"
  ( cd "$SIM" && env PATH="/nonexistent" PYTHONPATH="$SIM/src" \
      "$PY" -m pytest tests -q -p no:cacheprovider "$@" ) || fail=1
}

# GPU image: torch is present, the full tree runs.
run_one "GPU image selection (full tree)"
# CPU image: torch is absent by design, so the two torch modules are excluded exactly as the
# CPU build excludes them.
#
# WHAT THIS SECOND RUN DOES AND DOES NOT PROVE. It proves the CPU build's SELECTION collects
# and passes against the image's file set. It does NOT prove the selection survives torch being
# absent, because whatever machine you are running on almost certainly has torch: a module
# outside the ignore list that grows a torch import would pass here and fail in the image. Only
# the torch-free `test` job in CI, or the real CPU build, establishes that. Same one-directional
# blind spot ci_local.sh carries, and worth knowing before trusting a green line here.
run_one "CPU image selection (${#CPU_IGNORE[@]} module(s) excluded)" "${CPU_IGNORE[@]}"

if [ "$fail" -ne 0 ]; then
  echo
  echo "FAIL: the suite does not pass against the container's file set, so the image build"
  echo "will fail and the PUBLISHED IMAGE WILL SILENTLY STAY STALE. Fix the test to skip"
  echo "or degrade where the file it needs is absent, rather than removing it."
  exit 1
fi
echo
echo "OK: both image builds' test steps will pass."
