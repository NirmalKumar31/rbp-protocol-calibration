"""Every completion marker must land where its job's identity is allowed to write.

WHY THIS TEST EXISTS. Stage 11 ran assign, score and phyloP to completion, uploaded all
three result tables, and then failed with

    403 rbp-analysis@... does not have storage.objects.create access to
        .../buckets/rbp-repro-2026-derived/objects/variants-complete.json

Ninety minutes of correct work discarded by the last line of the stage, four times over. The
IAM condition was right: rbp-analysis is deliberately scoped to results/, variants/ and
driver/. The marker path was wrong -- the bucket root is outside every allowed prefix.

Stage 13 had the identical defect waiting, same identity and same root path, and would have
failed the same way after its own successful run.

A grep would not have caught this: both paths look perfectly ordinary in isolation. What makes
them wrong is a fact that lives in Terraform, so the check has to compare the two.
"""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]

# Which identity runs which stage script, from cloud/submit.sh's per-job service accounts.
SCRIPT_IDENTITY = {
    "scripts/cloud_variants.py": "rbp-analysis",
    "scripts/cloud_analysis.py": "rbp-analysis",
}

# The object prefixes each identity may create objects under, from cloud/terraform/iam.tf.
# Kept here as a literal on purpose: if someone widens the real condition, this test should
# fail and make them say so out loud rather than silently agreeing.
WRITABLE = {
    "rbp-analysis": ("results/", "variants/", "driver/"),
    "rbp-modal": ("runs/", "ckpt/", "variants/"),
    "rbp-train": ("runs/", "ckpt/", "rehearsal/", "results/"),
    "rbp-prep": ("processed/", "panel/", "manifest/", "interim/"),
}


def _marker(path):
    text = (ROOT / path).read_text()
    m = re.search(r'^MARKER\s*=\s*["\']([^"\']+)["\']', text, re.M)
    assert m, f"{path} defines no MARKER constant"
    return m.group(1)


@pytest.mark.parametrize("script,identity", sorted(SCRIPT_IDENTITY.items()))
def test_marker_is_under_a_writable_prefix(script, identity):
    marker = _marker(script)
    allowed = WRITABLE[identity]
    assert marker.startswith(allowed), (
        f"{script} writes its completion marker to {marker!r}, but {identity} may only "
        f"create objects under {allowed}. The stage will do all of its work and then 403 on "
        f"the final upload. Move the marker under one of those prefixes."
    )


@pytest.mark.parametrize("script,identity", sorted(SCRIPT_IDENTITY.items()))
def test_marker_is_not_at_the_bucket_root(script, identity):
    marker = _marker(script)
    assert "/" in marker, (
        f"{script} writes {marker!r} at the bucket root. No service account in this project "
        f"has an unconditional write on any bucket, by design, so a root-level object is "
        f"always a 403."
    )
