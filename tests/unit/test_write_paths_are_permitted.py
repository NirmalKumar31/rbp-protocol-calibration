"""Every GCS object a stage writes must sit under a prefix that stage's identity may write.

WHY THIS FILE EXISTS, and it is the third time this bug shipped. Buckets in this project are
not writable wholesale: each service account carries an IAM condition restricting
storage.objects.create to a few object prefixes. That is deliberate and it works. What keeps
going wrong is the other side of the contract -- code that writes somewhere outside its own
prefix list, which cannot fail until the moment of the upload.

The failure is maximally expensive by construction, because uploads happen at the END of a
stage:

  * stage 11 ran assign, score and phyloP over 66,010 variants, uploaded three tables, then
    403'd writing variants-complete.json at the bucket root. Ninety minutes, four times.
  * stage 13 had the identical defect queued behind it, same identity, same root path.
  * the fix for those introduced a THIRD, writing the phyloP cache to interim/ -- a prefix
    belonging to rbp-prep -- by analogy with the local filesystem path.
  * stage 12 would have cut and uploaded 164,835 variant windows and then 403'd writing its
    task manifest to manifest/, also rbp-prep's.

Three of those four were found by reading code against Terraform rather than by running a
job. This test makes that audit automatic, because the lesson from every previous round is
that a check which only fires at runtime fires too late.

WHAT THIS DOES NOT COVER. Only literal prefixes recoverable statically. A path assembled at
runtime from data cannot be checked here, so the convention is to keep the prefix literal in
the f-string -- `f"variants/{sub}/{name}.csv"` is checkable, `f"{dest}/{name}.csv"` is not.
"""

import ast
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]

# From cloud/terraform/iam.tf. Duplicated here on purpose: if someone widens a real condition,
# this test should fail and make them justify it rather than silently agreeing.
WRITABLE = {
    "rbp-analysis": ("results/", "variants/", "driver/"),
    "rbp-modal": ("runs/", "ckpt/", "variants/"),
    "rbp-train": ("runs/", "ckpt/", "rehearsal/", "results/"),
    "rbp-prep": ("processed/", "panel/", "manifest/", "interim/"),
}

# Which identities execute which file. Stage scripts run under their job's service account;
# the ones invoked on the driver VM inherit the driver's identity, which is rbp-analysis.
#
# variant_splicebert.py RUNS UNDER TWO IDENTITIES and that is the interesting case.
# `--what tables` runs on the driver as rbp-analysis; `--what task` runs inside the Modal
# container as rbp-modal, because modal_variants.py does not touch GCS itself -- it shells out
# to this same script with the Modal secret's key. So its writes must satisfy BOTH, and the
# intersection of (results/, variants/, driver/) and (runs/, ckpt/, variants/) is variants/
# alone. A path that is legal for one identity and not the other fails only on whichever half
# of the pipeline runs second, which is the worst possible time to find out.
RUNS_AS = {
    "scripts/cloud_variants.py": ("rbp-analysis",),
    "scripts/cloud_analysis.py": ("rbp-analysis",),
    "scripts/variant_splicebert.py": ("rbp-analysis", "rbp-modal"),
    "cloud/modal/modal_variants.py": ("rbp-modal",),
}


def _allowed(identities):
    """Prefixes permitted to EVERY identity that runs the file."""
    sets = [set(WRITABLE[i]) for i in identities]
    return tuple(sorted(set.intersection(*sets)))

# Calls that create or overwrite an object. Reads are unrestricted for these identities.
WRITE_CALLS = ("upload_from_filename", "upload_from_string", "upload_from_file")


def _string_assignments(tree):
    """NAME -> leading literal text, for every string assignment anywhere in the module.

    Function-local assignments are included, not just module constants. The real keys are
    built as locals -- `key = f"variants/tables/{cell}_{prot}.tsv.gz"` then
    `bucket.blob(key).upload_from_string(...)` -- so a module-level-only scan reported them as
    unresolvable and skipped them. Two of the four files were being waved through.

    Name collisions across scopes are possible in principle; a conflicting reassignment would
    make this stricter rather than laxer, since any one of the recorded values failing the
    prefix check fails the test.
    """
    out = {}
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        lit = None
        if isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
            lit = node.value.value
        elif isinstance(node.value, ast.JoinedStr) and node.value.values \
                and isinstance(node.value.values[0], ast.Constant):
            lit = node.value.values[0].value
        if lit is None:
            continue
        for t in node.targets:
            if isinstance(t, ast.Name):
                out.setdefault(t.id, []).append(lit)
    return out


def _static_prefixes(node, names):
    """Every literal leading text a blob key could have. Empty if none can be known."""
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return [node.value]
    if isinstance(node, ast.Name):
        return list(names.get(node.id, []))
    if isinstance(node, ast.JoinedStr):        # f-string: take the leading literal run
        if node.values and isinstance(node.values[0], ast.Constant):
            return [node.values[0].value]
        return []
    return []


def _written_keys(path):
    """(line, key) for every blob key a write call is invoked on. key None if unresolvable."""
    tree = ast.parse((ROOT / path).read_text())
    names = _string_assignments(tree)
    found = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        f = node.func
        if not isinstance(f, ast.Attribute) or f.attr not in WRITE_CALLS:
            continue
        # .blob(KEY).upload_from_*(...)  ->  walk back to the blob() call
        inner = f.value
        while isinstance(inner, ast.Call) and isinstance(inner.func, ast.Attribute) \
                and inner.func.attr not in ("blob",):
            inner = inner.func.value
        if isinstance(inner, ast.Call) and isinstance(inner.func, ast.Attribute) \
                and inner.func.attr == "blob" and inner.args:
            keys = _static_prefixes(inner.args[0], names)
            found += [(node.lineno, k) for k in keys] or [(node.lineno, None)]
    return found


@pytest.mark.parametrize("path,identities", sorted(RUNS_AS.items()))
def test_every_static_write_is_under_a_permitted_prefix(path, identities):
    if not (ROOT / path).exists():
        pytest.skip(f"{path} absent in this checkout")
    allowed = _allowed(identities)
    bad = []
    for lineno, key in _written_keys(path):
        if key is None:
            continue                            # runtime-assembled, see module docstring
        if not key.startswith(allowed):
            bad.append(f"{path}:{lineno} writes {key!r}")
    assert not bad, (
        f"{'+'.join(identities)} may only create objects under {allowed}, but:\n  "
        + "\n  ".join(bad)
        + "\nThese fail with 403 storage.objects.create at the END of the stage, after all "
          "the work is done. Move the object under a permitted prefix, or change which "
          "identity runs the stage."
    )


@pytest.mark.parametrize("path,identities", sorted(RUNS_AS.items()))
def test_no_write_lands_at_the_bucket_root(path, identities):
    if not (ROOT / path).exists():
        pytest.skip(f"{path} absent in this checkout")
    root_writes = [f"{path}:{ln} writes {k!r}"
                   for ln, k in _written_keys(path) if k and "/" not in k]
    assert not root_writes, (
        "No identity in this project has an unconditional write on any bucket, by design, so "
        "a root-level object is always a 403:\n  " + "\n  ".join(root_writes))


def test_the_iam_table_here_matches_terraform():
    """If Terraform's conditions drift from this table, every check above is worthless.

    Skipped inside the container, which ships src/, scripts/, config/ and tests/ but not
    cloud/ -- there is no reason for a worker to carry Terraform. The image's own test gate
    caught this as a FileNotFoundError, which is the gate doing its job.
    """
    iam = ROOT / "cloud" / "terraform" / "iam.tf"
    if not iam.exists():
        pytest.skip("cloud/terraform not present (running inside the container image)")
    tf = iam.read_text()
    for sa, prefixes in WRITABLE.items():
        for p in prefixes:
            assert f'/objects/{p}"' in tf, (
                f"{sa}: prefix {p!r} is listed here but no matching startsWith() appears in "
                f"iam.tf. Either Terraform changed or this table is stale.")
