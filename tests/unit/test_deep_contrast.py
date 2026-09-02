"""The R1g transplant must be the same arithmetic as R1b's, and must not flatter itself.

scale_check.py and deep_model_contrast.py both split a contrast into the part AUROC
compression explains and the part it does not. They are separate files, so the second can
drift from the first without anything failing -- and the whole force of R1g is that the deep
models were put through the SAME decomposition as the k-mer, not a friendlier one.

These tests pin the properties that make the split honest rather than the numbers, which
live in golden.yaml.
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def dmc():
    return _load("dmc", "scripts/deep_model_contrast.py")


@pytest.fixture(scope="module")
def sc():
    return _load("sc", "scripts/scale_check.py")


def test_dprime_and_auroc_are_inverses(dmc):
    a = np.array([0.51, 0.60, 0.75, 0.90, 0.99])
    assert np.allclose(dmc.auroc(dmc.dprime(a)), a)


def test_transplant_matches_scale_check_arithmetic(dmc, sc):
    """The d' transform must be identical in both files, not merely similar."""
    a = np.array([0.55, 0.62, 0.78, 0.91])
    assert np.allclose(dmc.dprime(a), sc.dprime(a))
    assert np.allclose(dmc.auroc(dmc.dprime(a)), sc.auroc(sc.dprime(a)))


def _frame(comp_gc, comp_dn, full_gc, full_dn):
    d = pd.DataFrame({"comp_gc": comp_gc, "comp_dn": comp_dn,
                      "m_full_gc": full_gc, "m_full_dn": full_dn})
    d["m_gain_gc"] = d.comp_gc.rsub(d.m_full_gc)
    d["m_gain_dn"] = d.comp_dn.rsub(d.m_full_dn)
    return d


def test_no_protocol_effect_when_only_the_baseline_differs(dmc):
    """THE NULL THIS SPLIT EXISTS TO REPRESENT.

    Give both arms the same d' increment on different baselines. The whole contrast is then
    compression by construction, and the protocol effect must be zero -- in BOTH transplant
    directions. A split that returned something positive here would manufacture R1g out of
    the ceiling alone.
    """
    comp_gc = np.array([0.78, 0.82, 0.70])
    comp_dn = np.array([0.62, 0.65, 0.58])
    dd = 0.45                                   # one shared increment, d' units
    full_gc = dmc.auroc(dmc.dprime(comp_gc) + dd)
    full_dn = dmc.auroc(dmc.dprime(comp_dn) + dd)
    q = dmc.transplant(_frame(comp_gc, comp_dn, full_gc, full_dn), "m")
    assert abs(q["contrast_protocol"]) < 1e-12
    assert abs(q["contrast_protocol_reverse"]) < 1e-12
    assert q["contrast_auroc"] == pytest.approx(q["contrast_scale_only"], abs=1e-12)


def test_protocol_effect_appears_only_with_a_real_increment_gap(dmc):
    comp_gc = np.array([0.78, 0.82, 0.70])
    comp_dn = np.array([0.62, 0.65, 0.58])
    full_gc = dmc.auroc(dmc.dprime(comp_gc) + 0.30)
    full_dn = dmc.auroc(dmc.dprime(comp_dn) + 0.60)   # dinuc genuinely gains more
    q = dmc.transplant(_frame(comp_gc, comp_dn, full_gc, full_dn), "m")
    assert q["contrast_protocol"] > 0
    assert q["contrast_protocol_reverse"] > 0


def test_transplant_family_has_four_members_and_is_reported_as_a_range(dmc):
    """Reporting one member would be question-begging; the script must emit all four."""
    comp_gc = np.array([0.78, 0.82])
    comp_dn = np.array([0.62, 0.65])
    full_gc = dmc.auroc(dmc.dprime(comp_gc) + 0.30)
    full_dn = dmc.auroc(dmc.dprime(comp_dn) + 0.60)
    q = dmc.transplant(_frame(comp_gc, comp_dn, full_gc, full_dn), "m")
    fam = ("contrast_protocol", "contrast_protocol_reverse",
           "contrast_protocol_logit", "contrast_protocol_logit_reverse")
    assert all(k in q for k in fam)
    assert len({round(q[k], 12) for k in fam}) > 1     # they genuinely differ


def test_min_coverage_floor_would_catch_a_missing_fold(dmc):
    """A fold that failed to upload costs ~20% of the rows; the floor must reject that."""
    assert dmc.MIN_COVERAGE > 0.8
    assert dmc.MIN_COVERAGE <= 0.999


def test_kmer_is_on_the_ladder(dmc):
    """R1g's comparison is only meaningful if the k-mer is measured the same way."""
    assert dmc.MODELS[0] == "kmer"
