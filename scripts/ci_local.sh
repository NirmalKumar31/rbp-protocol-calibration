#!/usr/bin/env bash
# Run EXACTLY what .github/workflows/ci.yml runs, locally, before pushing.
#
# WHY THIS EXISTS. CI failed four times on pushes that passed locally, and every time the
# reason was that the two ran different commands. The workflow excludes two test modules
# (they import torch, which is not in requirements-cpu.txt), which changes the collected
# count and therefore what tests/unit/test_suite_size.py asserts; it also runs
# cloud/package_repo.sh, which nothing local exercised. So "the suite passes" locally was
# never a prediction about CI, and the gap was rediscovered by email each time.
#
# The rule this encodes: if a check runs in CI it must be runnable in one command here, and
# this script is the single place the invocation is written down. When the workflow changes,
# change it here in the same commit.
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

# Mirrors the workflow's env block. GOOGLE_CLOUD_PROJECT is resolved at import by
# rbp.utils.cloud and raises with no default on purpose, so CI supplies a fake one.
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

# THIS IS NOW IN THE WORKFLOW TOO. It used to carry the note "not in the workflow, because it
# needs the committed tables and takes longer than CI should", which was wrong on both counts:
# the tables are committed, so a clean clone has them, and the run takes under a minute. The
# effect of leaving it out was that README's headline instruction was gated nowhere.
step "manuscript numbers trace to a table (BEFORE the verifier, which reads its output)"
"$PY" scripts/audit_manuscript.py | tail -3 || fail "audit_manuscript.py"

step "verifier"
"$PY" scripts/verify.py --local results/tables | tail -2 || fail "verify.py"

step "release documents are consistent with the artefacts"
"$PY" scripts/release_consistency.py || fail "release_consistency.py"

step "every committed table has a producing script"
"$PY" scripts/provenance.py --check || fail "provenance.py"

step "regenerated artefacts match what is committed"
git diff --exit-code -- results/tables/ >/dev/null || fail "a generated table changed; commit it"
echo "  clean"

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
