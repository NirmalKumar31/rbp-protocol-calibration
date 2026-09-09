#!/usr/bin/env bash
# Run the local release checks in dependency order before pushing.
#
#   ./scripts/preflight.sh            # check only; fails like CI would
#   ./scripts/preflight.sh --fix      # also sync derived counts and refresh manifests
#
# audit_manuscript.py must precede verify.py, and the two manifest generators must run last.
# In --fix mode, release metadata is synchronised before the checks run.
set -uo pipefail
cd "$(dirname "$0")/.."
PY="${PY:-python3}"
# Preserve any caller-supplied Python path after the project source directory.
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
export GOOGLE_CLOUD_PROJECT="${GOOGLE_CLOUD_PROJECT:-ci-no-such-project}"

FIX=0
[ "${1:-}" = "--fix" ] && FIX=1

# Run a declared subset when torch is unavailable; CI covers the complete environment.
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

# 1. Manuscript audit (verify.py reads its output).
step "manuscript numbers trace to a table" "$PY" scripts/audit_manuscript.py

# 2. Derived release metadata.
if [ "$FIX" = 1 ]; then
  step "release documents (syncing counts)" "$PY" scripts/release_consistency.py --fix
fi
# A full test census requires torch.
if [ "$HAVE_TORCH" = 1 ]; then
  step "release documents are consistent" "$PY" scripts/release_consistency.py --require-all
else
  step "release documents are consistent (partial: no test census)" \
       "$PY" scripts/release_consistency.py
fi

# In --fix mode, refresh once before tests and once after all generators.
if [ "$FIX" = 1 ]; then
  step "refresh manifests (early pass, so the suite sees them fresh)" \
       env PY="$PY" ./scripts/refresh_manifests.sh
fi

# 3. Published values.
step "published values verify offline" "$PY" scripts/verify.py --local results/tables

# 4. Documented offline entry points.
step "every --from-cache entry point reproduces" "$PY" scripts/cache_idempotence.py

# 5. Offline manifests.
step "raw input manifest" "$PY" scripts/raw_inputs.py --check
step "history secret scan" "$PY" scripts/history_scan.py --check

# 6. Code.
if [ "$HAVE_TORCH" = 1 ]; then
  step "unit suite" "$PY" -m pytest tests
else
  step "unit suite, minus the two torch modules" "$PY" -m pytest tests \
       --ignore=tests/unit/test_models.py --ignore=tests/unit/test_train_folds.py
  skip "the two torch test modules" \
       "no torch. Install 'pip install -e .[dev,neural] -c constraints.txt' to cover them"
fi
step "ruff" "$PY" -m ruff check .
# Fall back to find in exported archives, and fail if no shell scripts are found.
step "shell syntax" bash -c '
  files=$(git ls-files "*.sh" 2>/dev/null || find . -name "*.sh" -not -path "./.git/*")
  [ -n "$files" ] || { echo "no shell scripts found; this step checked NOTHING" >&2; exit 1; }
  for f in $files; do bash -n "$f" || exit 1; done'

# 6b. Test the reduced file set copied into the container image.
if [ "$HAVE_TORCH" = 1 ]; then
  step "the suite passes against the container's file set" \
       env PY="$PY" bash scripts/check_image_tree.sh
else
  skip "the container file-set check" \
       "it exercises the GPU image selection, which needs torch"
fi

# 6c. A hyphen ending a source line splits a word, and the render can hide it. Needs no
# toolchain, so it is unconditional.
step "no word is split across a line by a hyphen" "$PY" scripts/tex_line_breaks.py

# 6c. Check figure-text geometry when Poppler is available.
if command -v pdftotext >/dev/null 2>&1; then
  step "no figure draws text on top of text" "$PY" scripts/figure_overlap.py
else
  skip "the figure text-overlap check" \
       "no pdftotext (poppler-utils) here; the CI manuscript job installs it"
fi

# 7. Build in a temporary copy and compare the tracked PDFs with the result.
if command -v pdflatex >/dev/null 2>&1; then
  step "tracked PDFs match a clean build (builds in a temp copy)" "$PY" scripts/pdf_freshness.py
else
  skip "tracked PDFs match a clean build" \
       "no pdflatex here, so the CI manuscript job is the ONLY thing checking it"
fi

# 8. Refresh or validate manifests last.
if [ "$FIX" = 1 ]; then
  step "refresh manifests (two passes, LAST)" env PY="$PY" ./scripts/refresh_manifests.sh
else
  step "the column dictionary matches the tables" "$PY" scripts/column_dictionary.py --check
  step "every committed table has a producing script" "$PY" scripts/provenance.py --check
fi

printf '\n'
if [ "$fail" = 0 ] && [ "$skipped" = 1 ]; then
  printf 'PREFLIGHT PARTIAL. Everything runnable here passed, but the steps marked SKIPPED\n'
  printf 'above did NOT run, so this is NOT the full release gate. For that: install\n'
  printf 'pip install -e ".[dev,neural]" -c constraints.txt and a TeX distribution, or read\n'
  printf 'the GitHub Actions run for this commit, whose jobs cover both.\n'
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
