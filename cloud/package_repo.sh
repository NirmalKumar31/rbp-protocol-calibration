#!/usr/bin/env bash
# Package the repo for the driver VM, and PROVE the package is complete before uploading.
#
# WHY THIS IS A SCRIPT AND NOT A TAR COMMAND. It was a tar command, typed inline:
#
#   tar czf repo.tgz --exclude='.git' --exclude='data' --exclude='results' cloud scripts src ...
#
# `--exclude='data'` matches a path COMPONENT at any depth, so it silently dropped
# src/rbp/data/ -- a source package, not a data directory. The driver unpacked a repo whose
# `rbp.data` module did not exist, Modal uploaded that src tree to 188 containers, and every
# one died with ModuleNotFoundError: No module named 'rbp.data' before doing any work.
#
# This is the SAME BUG the .gitignore had, in a different tool. There, `data/` matched
# src/rbp/data/ and shipped an image without a source package; the fix was to anchor it as
# `/data/`. Knowing that story did not stop me typing the unanchored version into tar.
#
# The excludes were pointless as well as harmful: the tar list names explicit directories and
# none of them is a top-level data/ or results/. So they are gone, and what replaces them is
# a check that every package present in src/rbp/ is present in the archive.

set -uo pipefail
cd "$(dirname "$0")/.." || exit 1

OUT="${1:-/tmp/repo.tgz}"
INCLUDE=(cloud scripts src config tests docker run.sh pyproject.toml)

# Anchored, so they cannot match a nested path component. --exclude=./x with tar's leading-./
# form, plus the vcs/cache patterns which are safe because they are unambiguous names.
tar czf "$OUT" \
    --exclude='__pycache__' --exclude='*.pyc' --exclude='.terraform' --exclude='tfplan*' \
    "${INCLUDE[@]}" || { echo "tar failed" >&2; exit 1; }

# --- the verification, which is the entire point of the file -----------------------------
fail=0

# Every package directory under src/rbp/ must appear in the archive. Compared against the
# working tree rather than a hardcoded list, so a new package is covered the day it is added.
for d in src/rbp/*/; do
  name=$(basename "$d")
  [ "$name" = "__pycache__" ] && continue
  if ! tar tzf "$OUT" | grep -q "^src/rbp/${name}/"; then
    echo "MISSING PACKAGE: src/rbp/${name}/ is in the working tree but not in $OUT" >&2
    fail=1
  fi
done

# The files the driver actually invokes. A tarball that unpacks but cannot run is no better
# than a missing one.
for f in cloud/submit.sh cloud/driver_startup.sh scripts/cloud_variants.py \
         scripts/cloud_analysis.py scripts/variant_splicebert.py \
         cloud/modal/modal_variants.py cloud/modal/modal_sweep.py; do
  tar tzf "$OUT" | grep -qx "$f" || { echo "MISSING FILE: $f" >&2; fail=1; }
done

# Importability is the property that actually matters, so test it rather than infer it.
# Unpack into a scratch dir and import every module the driver and Modal rely on.
TMP=$(mktemp -d) || exit 1
trap 'rm -rf "$TMP"' EXIT
tar xzf "$OUT" -C "$TMP"
PY="${PY:-python3}"
GOOGLE_CLOUD_PROJECT="${GOOGLE_CLOUD_PROJECT:-packaging-check}" \
PYTHONPATH="$TMP/src" "$PY" - <<'EOF' || fail=1
import importlib, sys
mods = ["rbp", "rbp.data", "rbp.data.windows", "rbp.variants", "rbp.variants.assign",
        "rbp.variants.phylop", "rbp.variants.clinvar", "rbp.eval", "rbp.utils.cloud"]
bad = []
for m in mods:
    try:
        importlib.import_module(m)
    except Exception as e:
        bad.append(f"{m}: {type(e).__name__}: {e}")
if bad:
    print("IMPORT FAILURES IN THE PACKAGED TREE:", file=sys.stderr)
    for b in bad:
        print("  " + b, file=sys.stderr)
    sys.exit(1)
print(f"  all {len(mods)} modules import from the packaged tree")
EOF

[ "$fail" -eq 0 ] || { echo "package verification FAILED; not uploading" >&2; exit 1; }
echo "packaged $OUT ($(du -h "$OUT" | cut -f1)), $(tar tzf "$OUT" | wc -l | tr -d ' ') entries, verified"
