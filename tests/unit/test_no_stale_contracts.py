"""Present-tense claims in prose and comments must agree with the released evidence.

WHY. An external audit found six places where the repository asserts, in the present tense,
something a committed table or another file in the same repository contradicts. None of them
changed a number. All of them would be found by a referee, and each one costs more credibility
than the defect it describes, because a reader cannot tell a stale comment from a live one.

The shape recurs: a disclosure is written carefully, the defect is then fixed, and the
disclosure stays. `deep_model_contrast.py` disclosed fold leakage in 20 datasets and ended
with "do not restore the same-folds sentence without rerunning the sweep". The sweep was rerun
in 2276eea, `fold_integrity.csv` went from 20 to 0, and the disclosure sat there for three
release candidates telling the reader that a fifth of one arm was contaminated.

Each check below pins one claim to the artefact that decides it, so the claim cannot outlive
the state it describes. These are string checks, which is crude, but the failure they catch is
also crude and nothing subtler was catching it.
"""

import csv
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SECTIONS = ROOT / "manuscript" / "sections"


def _read(rel):
    p = ROOT / rel
    if not p.exists():
        pytest.skip(f"{rel} not in this checkout")
    return p.read_text()


def test_the_discussion_does_not_deny_an_experiment_the_results_report():
    """P0.1. Both statements were in the same PDF, 850 lines apart."""
    results = _read("manuscript/sections/results.tex")
    disc = _read("manuscript/sections/discussion.tex")
    if "\\label{sec:transport}" not in results:
        pytest.skip("no transport section to be consistent with")
    assert "performs that transfer and we do not" not in disc, (
        "results.tex reports the full train-arm by evaluation-arm experiment under "
        "\\label{sec:transport}, and discussion.tex says the paper did not perform it")
    assert "\\ref{sec:transport}" in disc, (
        "the Discussion's scope paragraph must point at the transport section rather than "
        "describe the design as if it were absent")


def test_the_model_was_refit_per_arm_so_nothing_says_the_model_was_held_fixed():
    """P0.2. The model class was fixed. A separate model was fitted in every arm."""
    bad = re.compile(r"[Hh]olding the model,|\bModel, peak set\b")
    for rel in ("manuscript/sections/introduction.tex", "manuscript/sections/discussion.tex",
                "CITATION.cff"):
        t = _read(rel)
        assert not bad.search(t), (
            f"{rel} says the MODEL was held fixed across arms. It was refit in each arm; what "
            "was held fixed is the model class and its hyperparameters")


def test_the_gc_and_dinucleotide_positives_are_not_called_identical():
    """P0.2. Only 10 of 94 datasets have identical positive sets."""
    t = _read("manuscript/sections/results.tex")
    assert "identical positives" not in t, (
        "the GC and dinucleotide arms retain near-identical, not identical, positives: each "
        "arm keeps only the positives its own matching could pair")
    overlap = ROOT / "results" / "tables" / "positive_set_overlap.csv"
    if not overlap.exists():
        return
    j = [float(r["jaccard"]) for r in csv.DictReader(overlap.open())]
    assert min(j) < 1.0, (
        "every dataset now has identical GC and dinucleotide positives, so the caption may "
        "say identical again and this test should be retired")
    assert sum(x == 1.0 for x in j) < len(j), "same"
    # The caption quotes these three. Pin them so a regenerated table cannot silently
    # contradict the text.
    assert abs(sorted(j)[len(j) // 2] - 0.9972) < 5e-4, "median Jaccard moved"
    assert abs(min(j) - 0.9164) < 5e-4, "minimum Jaccard moved"
    assert sum(x == 1.0 for x in j) == 10, "the count of identical datasets moved"


def test_the_fold_leakage_disclosure_matches_the_fold_integrity_table():
    """P1.13. The disclosure outlived the defect by three release candidates."""
    fi = ROOT / "results" / "tables" / "fold_integrity.csv"
    if not fi.exists():
        pytest.skip("fold_integrity.csv not in this checkout")
    rows = {r["check"]: float(r["value"]) for r in csv.DictReader(fi.open())}
    dirty = rows.get("datasets NOT chromosome-grouped, dn arm")
    src = _read("scripts/deep_model_contrast.py")
    present = re.search(r"For \*\*20 of the 94 dinucleotide-arm\s*\n?datasets\*\* the\s*\n"
                        r"committed", src)
    if dirty == 0:
        assert present is None, (
            "fold_integrity.csv reports 0 non-chromosome-grouped dn-arm datasets, and "
            "deep_model_contrast.py still asserts in the present tense that 20 of them carry "
            "scores from a stratified random partition")
    else:
        assert present is not None, (
            f"fold_integrity.csv reports {dirty} non-chromosome-grouped dn-arm datasets and "
            "deep_model_contrast.py no longer discloses it")


def test_the_baseline_is_not_called_the_whole_story():
    """P1.13. The same script prints that a protocol-specific residual exists."""
    src = _read("scripts/protocol_or_baseline.py")
    assert "A protocol-specific residual DOES exist" in src, (
        "this test pins the docstring to that finding; if the finding moved, move the test")
    for phrase in ("it is the whole story", "carries essentially no information beyond it"):
        assert phrase not in src, (
            f"protocol_or_baseline.py's header says {phrase!r} while the script itself reports "
            "a residual whose interval excludes zero under two designs")


def test_the_killswitch_does_not_promise_that_nothing_is_lost():
    """P1.11. Terraform promised what the function's own docstring retracts."""
    tf = _read("cloud/terraform/killswitch.tf")
    assert "never lost data" not in tf, (
        "killswitch.tf guarantees no data loss on billing detach. Google documents that "
        "disabling billing may delete resources non-recoverably, and cloud/killswitch/main.py "
        "says so. An unsafe guarantee is what persuades an operator to skip the backup")


def test_the_venue_is_not_described_as_arxiv():
    """P1.12. The target has been bioRxiv since the venue was chosen."""
    assert "arXiv PDF" not in _read("manuscript/build.sh")


def test_exit_zero_is_not_described_as_raw_reproduction():
    """P1.7/P0.4. Most released tables are verified, not reconstructed."""
    t = _read("docs/REPRODUCE.md")
    assert "Exit 0 means the science reproduced." not in t, (
        "run.sh regression-verifies most tables against committed evidence rather than "
        "rebuilding them from raw inputs; PROVENANCE.csv says which is which")


def test_the_image_gate_resolves_an_absolute_interpreter():
    """P0.4. `env PATH=/nonexistent python3` can never find python3.

    run.sh gates the paid image build on this script, so a gate that cannot start is a rebuild
    path that is documented and does not run. It failed for everyone, every time, with
    `env: python3: No such file or directory`, and nothing noticed because nothing ran it.
    """
    t = _read("scripts/check_image_tree.sh")
    assert 'PY="${PY:-python3}"' not in t, (
        "PY must be resolved to an absolute path with command -v BEFORE PATH is emptied")
    assert "command -v" in t, "the interpreter is not being resolved at all"


def test_the_image_gate_runs_both_build_selections():
    """P0.4. One simulation cannot stand in for two different images."""
    t = _read("scripts/check_image_tree.sh")
    assert "cloudbuild.cpu.yaml" in t, (
        "the CPU build's --ignore list must be read from the build file, not copied here, "
        "or the simulation drifts from the thing it simulates")
    cpu = _read("docker/cloudbuild.cpu.yaml")
    ignores = re.findall(r"--ignore=tests/[A-Za-z0-9_/.]+", cpu)
    assert ignores, (
        "the parser in check_image_tree.sh finds the CPU ignore list with this pattern; if "
        "cloudbuild.cpu.yaml stops matching it, the gate fails closed but for the wrong reason")


def test_the_neural_stack_is_an_extra_not_a_base_dependency():
    """P0.4. `pip install -e . -c constraints.txt` used to install torch regardless."""
    tomllib = pytest.importorskip("tomllib")
    d = tomllib.loads(_read("pyproject.toml"))
    base = " ".join(d["project"]["dependencies"])
    for pkg in ("torch", "transformers", "multimolecule", "peft"):
        assert pkg not in base, (
            f"{pkg} is a base dependency, so the offline verification install is not minimal. "
            "constraints.txt has always claimed it is")
    assert "neural" in d["project"]["optional-dependencies"]
    rp = d["project"]["requires-python"]
    assert "<" in rp, f"requires-python {rp!r} has no upper bound"
