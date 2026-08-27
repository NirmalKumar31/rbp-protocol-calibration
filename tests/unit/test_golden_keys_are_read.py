"""Every key in golden.yaml must actually be read by something. Three strikes.

WHY THIS TEST EXISTS, and it is the most embarrassing entry in the project.

`config/golden.yaml` is the file that makes this pipeline's claims checkable. A key sitting in
it unread is worse than no key at all: it reads as a guarantee, it survives review because
reviewers read the config and not the reader, and it makes the verifier's pass count a lie by
omission. This happened three times.

  1. `integrity.min_tests_passing` sat unread while the suite grew from 480 to 576, so the
     floor it enforced was four hundred tests stale.
  2. `r1_headline_is_gc_share_only` was added to FORBID a headline. Nothing read it, and the
     forbidden headline was duly promoted two weeks later -- by the person who wrote the key.
  3. The fix for (1) and (2) wired up 26 keys and, in the same commit, added a 9-key
     `strand_audit` block that nothing read, plus a line in `docs/59` asserting it was gated.

Three occurrences of one bug, each found by a human squinting at a grep. So the defence is
mechanical: this test parses golden.yaml, walks every leaf key, and fails if the name appears
nowhere in the code that consumes it. Inspection does not scale; a failing build does.

WHAT COUNTS AS "READ". The key's name appearing as a string in `scripts/verify.py` or in a
test. That is a weak check -- it cannot tell a real `spec["x"]` from the word x in a comment --
but it is strong enough for the failure mode that actually occurred, which was a key existing
in no file but the config. Anything stronger (AST-resolving `spec[...]` subscripts) breaks on
the loop-driven checks where the key name is built from a tuple, and a test that forces the
production code to be written in a less readable way to satisfy the test is the wrong trade.

DOCUMENTED EXEMPTIONS ONLY. `meta.*` keys are provenance for a human reader and are asserted
nowhere by design. Every other exemption needs a reason written next to it, which is the point:
adding a key to this set is a visible act in a diff, and forgetting to read one is not.
"""

from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
GOLDEN = ROOT / "config" / "golden.yaml"

# Files that legitimately consume golden keys.
READERS = ("scripts/verify.py",
           "tests/unit/test_suite_size.py",
           "tests/unit/test_golden_keys_are_read.py")

# Keys that are documentation for a human, not assertions. Each needs a stated reason.
EXEMPT = {
    "meta.reference_run": "provenance string, printed not asserted",
    "meta.established": "provenance string, printed not asserted",
    "meta.study_panel_n": "cross-reference for a reader; the count is asserted from the "
                          "tables themselves in test_panel_counts.py",
    "meta.r1_panel_n": "same",
    "meta.stage7_tasks": "same; the task count is derived from the manifest at runtime",
    "r2_four_models.required_order": "a list, consumed positionally by verify_r2's ordering "
                                     "check rather than by name",
    "donor_overlap_RETRACTED.see": "a pointer to the block that replaced it; the retraction "
                                   "is the content",
}


def _leaf_keys(obj, prefix=""):
    """Every scalar-or-{value,tol} leaf, as a dotted name."""
    if isinstance(obj, dict):
        # a {value:, tol:} pair is the leaf, not its two members
        if set(obj) <= {"value", "tol"} and "value" in obj:
            yield prefix
            return
        for k, v in obj.items():
            yield from _leaf_keys(v, f"{prefix}.{k}" if prefix else k)
    else:
        yield prefix


@pytest.fixture(scope="module")
def golden():
    return yaml.safe_load(GOLDEN.read_text())


@pytest.fixture(scope="module")
def reader_text():
    out = []
    for r in READERS:
        p = ROOT / r
        if p.exists():
            out.append(p.read_text())
    assert out, f"none of the reader files exist: {READERS}"
    return "\n".join(out)


def test_golden_file_parses(golden):
    assert isinstance(golden, dict) and golden, "golden.yaml is empty or not a mapping"


def test_every_key_is_read(golden, reader_text):
    """The whole point of the file. A key nobody reads is a guarantee nobody checked."""
    unread = []
    for dotted in _leaf_keys(golden):
        if dotted in EXEMPT:
            continue
        leaf = dotted.split(".")[-1]
        if leaf not in reader_text:
            unread.append(dotted)
    assert not unread, (
        f"{len(unread)} golden key(s) are read by nothing. Either wire them into "
        f"scripts/verify.py or add them to EXEMPT with a reason:\n  "
        + "\n  ".join(unread))


def test_exemptions_still_exist(golden):
    """An exemption for a key that has been deleted is stale permission. Fail on it."""
    present = set(_leaf_keys(golden))
    stale = [k for k in EXEMPT if k not in present]
    assert not stale, f"EXEMPT names keys that no longer exist in golden.yaml: {stale}"


def test_exemptions_have_reasons():
    assert all(v.strip() for v in EXEMPT.values()), "every exemption needs a stated reason"
