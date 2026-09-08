#!/usr/bin/env bash
# Mirror the GitHub Actions checks locally before pushing. Keep this file aligned with
# .github/workflows/ci.yml.
#
#   ./scripts/ci_local.sh
#
# Exits non-zero on the first failure, like CI. Add --fast to skip the package check, which
# is the slow step, when iterating on tests only.
set -uo pipefail
cd "$(dirname "$0")/.."

PY="${PY:-python3}"
FAST=0
[ "${1:-}" = "--fast" ] && FAST=1

step() { printf '\n=== %s\n' "$*"; }
fail() { printf 'CI-LOCAL FAILED: %s\n' "$*" >&2; exit 1; }

# Mirror the workflow environment without contacting GCP.
export GOOGLE_CLOUD_PROJECT="ci-no-such-project"
export PYTHONPATH="src"

step "unit tests (the workflow's exclusions, which change the collected count)"
"$PY" -m pytest tests/unit -q \
  --ignore=tests/unit/test_models.py \
  --ignore=tests/unit/test_train_folds.py \
  || fail "unit tests"

step "golden values checksum"
"$PY" - <<'EOF' || exit 1
import hashlib, pathlib
h = hashlib.sha256(pathlib.Path("config/golden.yaml").read_bytes()).hexdigest()
print("golden.yaml sha256:", h)
EOF

if [ "$FAST" = "0" ]; then
  step "repo package verifies"
  PY="$PY" ./cloud/package_repo.sh /tmp/ci_local_repo.tgz >/dev/null \
    || fail "cloud/package_repo.sh"
  echo "  package verified"
fi

step "manuscript numbers trace to a table (BEFORE the verifier, which reads its output)"
"$PY" scripts/audit_manuscript.py | tail -3 || fail "audit_manuscript.py"

step "verifier"
"$PY" scripts/verify.py --local results/tables | tail -2 || fail "verify.py"

step "release documents are consistent with the artefacts (full-suite job: --require-all)"
"$PY" scripts/release_consistency.py --require-all || fail "release_consistency.py"

# This command uses the current environment. Check Actions separately for the torch-free job.
echo "  NOTE: run with torch present, so this mirrors the full-suite job and not the"
echo "        torch-free 'test' job. Check the Actions log for that one."

step "the column dictionary is current"
"$PY" scripts/column_dictionary.py --check || fail "column_dictionary.py --check"

step "every committed table has a producing script"
"$PY" scripts/provenance.py --check || fail "provenance.py"

# Compare tracked PDFs with a clean build when TeX is available.
if command -v pdflatex >/dev/null 2>&1; then
  "$PY" scripts/pdf_freshness.py || fail "pdf_freshness.py"
else
  echo "  SKIPPED pdf_freshness.py: no pdflatex here, so the CI manuscript job is the only" >&2
  echo "  thing checking that the committed PDF matches the committed source" >&2
fi

step "regenerated artefacts match what is committed"
git diff --exit-code -- results/tables/ >/dev/null || fail "a generated table changed; commit it"
echo "  clean"

step "every --from-cache entry point reproduces its committed table"
"$PY" scripts/cache_idempotence.py || fail "cache_idempotence.py"

step "ruff"
if "$PY" -m ruff --version >/dev/null 2>&1; then
  "$PY" -m ruff check . || fail "ruff"
else
  echo "  SKIP: ruff not installed here (pip install ruff). CI runs it regardless."
fi

step "shell syntax"
for f in $(git ls-files '*.sh' 2>/dev/null || find . -name '*.sh' -not -path './.git/*'); do
  bash -n "$f" || fail "bash -n $f"
done
echo "  all shell scripts parse"

printf '\nALL CI-LOCAL CHECKS PASSED\n'
