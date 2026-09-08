#!/usr/bin/env bash
# Everything CI will check, in the order that works, before you push.
#
#   ./scripts/preflight.sh            # check only; fails like CI would
#   ./scripts/preflight.sh --fix      # also sync derived counts and refresh manifests
#
# Why this exists. 23 of the first 99 CI runs on this repository failed, and the failures are
# not random. Two families account for three quarters of them:
#
#   12 of 23   the unit suite, most often provenance: a table was regenerated and
#              refresh_manifests.sh was not run afterwards, or was run too early.
#    6 of 23   release-document consistency: an assertion was added, and the ~20 places that
#              state the assertion count went stale.
#
# Both are ORDERING and BOOKKEEPING failures, not defects. Neither is caught by running pytest,
# because pytest is not the thing that goes stale. They are caught by running the whole set in
# the right order, which is what CI does and what nobody does by hand at 1am.
#
# The two ordering rules that are not obvious, and that this script encodes so they cannot be
# violated by forgetting them:
#
#   1. audit_manuscript.py must run BEFORE verify.py, which reads manuscript_orphans.csv.
#      A stale orphan report once made the gate pass while a fresh audit found four.
#   2. refresh_manifests.sh must run LAST, after every other generator, and it is TWO PASSES
#      because the two manifest generators are circular. Anything regenerated afterwards
#      leaves the manifest stale again.
#
# --fix does the bookkeeping rather than reporting it: release_consistency --fix rewrites stale
# counts POSITIONALLY, at the offsets its own scan found, which is the only safe way to do it.
# Doing it with a regex over the digits is what rewrote fourteen unrelated scientific numbers on
# 2026-09-06 and 09-07.
set -uo pipefail
cd "$(dirname "$0")/.."
PY="${PY:-python3}"
# PREPEND, do not overwrite. Discarding a caller's PYTHONPATH is rude and it also made this
# script impossible to test against a simulated environment, which is how the torch-free path
# went unexercised.
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
export GOOGLE_CLOUD_PROJECT="${GOOGLE_CLOUD_PROJECT:-ci-no-such-project}"

FIX=0
[ "${1:-}" = "--fix" ] && FIX=1

# Is torch here? README's documented base install is `pip install -e . -c constraints.txt`,
# which deliberately omits torch, and this script then ran the whole test tree and the
# GPU-image test selection, both of which need it. So the documented commands could not
# succeed in the documented environment. An audit reproduced that in a clean archive.
#
# The answer is not to demand torch. It is to run what CAN run, say plainly what did not, and
# NOT print "PREFLIGHT CLEAN" over a subset, because a green word covering a partial run is the
# same defect as a gate nothing runs.
if "$PY" -c "import torch" >/dev/null 2>&1; then
  HAVE_TORCH=1
else
  HAVE_TORCH=0
fi
skipped=0
skip() {
  printf '\n=== %s\n    SKIPPED: %s\n' "$1" "$2" >&2
  skipped=1
}

fail=0
step() {
  local name="$1"; shift
  printf '\n=== %s\n' "$name"
  if "$@"; then
    printf '    OK\n'
  else
    printf '    FAILED: %s\n' "$name" >&2
    fail=1
  fi
}

# 1. the manuscript audit, first, because verify.py reads what it writes
step "manuscript numbers trace to a table" "$PY" scripts/audit_manuscript.py

# 2. the derived counts. With --fix this repairs them; without, it reports and fails.
if [ "$FIX" = 1 ]; then
  step "release documents (syncing counts)" "$PY" scripts/release_consistency.py --fix
fi
# --require-all demands every derived fact, and the test census cannot be derived without
# torch because the full suite will not collect. Without torch this reports rather than fails.
if [ "$HAVE_TORCH" = 1 ]; then
  step "release documents are consistent" "$PY" scripts/release_consistency.py --require-all
else
  step "release documents are consistent (partial: no test census)" \
       "$PY" scripts/release_consistency.py
fi

# TWICE, IN --fix, AND THE SECOND ONE IS STILL LAST. refresh_manifests.sh must run after every
# generator, which is why it is at the bottom. But tests/unit/test_provenance.py CHECKS those
# manifests and runs in the unit suite, well before it. So on any run where a table changed,
# --fix used to fail the provenance tests and pass on a second invocation, which trains you to
# run it twice and ignore the first result. Refreshing here as well means the suite sees fresh
# manifests; the pass at the bottom still catches anything the later steps regenerate.
if [ "$FIX" = 1 ]; then
  step "refresh manifests (early pass, so the suite sees them fresh)" \
       env PY="$PY" ./scripts/refresh_manifests.sh
fi

# 3. the published values
step "published values verify offline" "$PY" scripts/verify.py --local results/tables

# 4. every documented entry point returns its committed table
step "every --from-cache entry point reproduces" "$PY" scripts/cache_idempotence.py

# 5. the offline manifests a reader can check without the bucket
step "raw input manifest" "$PY" scripts/raw_inputs.py --check
step "history secret scan" "$PY" scripts/history_scan.py --check

# 6. code
if [ "$HAVE_TORCH" = 1 ]; then
  step "unit suite" "$PY" -m pytest tests
else
  step "unit suite, minus the two torch modules" "$PY" -m pytest tests \
       --ignore=tests/unit/test_models.py --ignore=tests/unit/test_train_folds.py
  skip "the two torch test modules" \
       "no torch. Install 'pip install -e .[neural] -c constraints.txt' to cover them"
fi
step "ruff" "$PY" -m ruff check .
step "shell syntax" bash -c 'for f in $(git ls-files "*.sh"); do bash -n "$f" || exit 1; done'

# 6b. THE CONTAINER'S FILE SET, which is smaller than this one and which nothing was running.
# docker/Dockerfile.cpu copies src, scripts, config, tests and pyproject.toml. Any test needing
# manuscript/ or results/ therefore FAILS inside the image, the image build fails, and the
# published image silently stays stale. check_image_tree.sh has existed to catch that since the
# image was unbuildable for weeks; it was in no pipeline, and by 2026-09-07 it was failing again
# on nine supplement tests added the day before. A gate nothing runs is not a gate.
if [ "$HAVE_TORCH" = 1 ]; then
  step "the suite passes against the container's file set" \
       env PY="$PY" bash scripts/check_image_tree.sh
else
  skip "the container file-set check" \
       "it exercises the GPU image selection, which needs torch"
fi

# 7. the manuscript, and whether the tracked PDF is the output of the tracked source
# pdf_freshness.py runs build.sh inside a TEMPORARY COPY, so every gate build.sh carries
# (undefined references, over- and underfull boxes) runs there and a failure propagates. It is
# therefore strictly stronger than an in-place build, and it does not rewrite the tracked PDFs.
# An in-place build here would dirty the tree on every check run, because pdflatex embeds a
# timestamp, and a "check" command that leaves the tree dirty is how timestamp churn gets
# committed by accident.
if command -v pdflatex >/dev/null 2>&1; then
  step "tracked PDFs match a clean build (builds in a temp copy)" "$PY" scripts/pdf_freshness.py
else
  printf '\n=== tracked PDFs match a clean build\n    SKIPPED: no pdflatex here.\n' >&2
  printf '    The CI manuscript job is then the ONLY thing checking this. That is a gap in\n' >&2
  printf '    this run, not an absence of one, and it is reported rather than passed over.\n' >&2
fi

# 8. LAST, and two passes. Anything regenerated after this leaves the manifest stale.
if [ "$FIX" = 1 ]; then
  step "refresh manifests (two passes, LAST)" env PY="$PY" ./scripts/refresh_manifests.sh
else
  step "the column dictionary matches the tables" "$PY" scripts/column_dictionary.py --check
  step "every committed table has a producing script" "$PY" scripts/provenance.py --check
fi

printf '\n'
if [ "$fail" = 0 ] && [ "$skipped" = 1 ]; then
  printf 'PREFLIGHT PARTIAL. Everything runnable here passed, but this environment has no\n'
  printf 'torch, so the steps marked SKIPPED above did NOT run and this is NOT the full\n'
  printf 'release gate. For that: pip install -e ".[neural]" -c constraints.txt, or read the\n'
  printf 'exact GitHub Actions run for this commit, whose full-suite job does have torch.\n'
elif [ "$fail" = 0 ]; then
  if [ "$FIX" = 1 ]; then
    printf 'PREFLIGHT CLEAN. Counts synced and manifests refreshed; review `git diff` and commit.\n'
  else
    printf 'PREFLIGHT CLEAN. This is what CI runs; if this passes, CI should too.\n'
  fi
else
  printf 'PREFLIGHT FAILED. Fix the steps above; pushing now will turn CI red.\n' >&2
  printf 'If the failures are stale counts or manifests, `./scripts/preflight.sh --fix`.\n' >&2
fi
exit "$fail"
