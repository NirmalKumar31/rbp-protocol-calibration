"""Weight initialisation must be seeded, and the training entry points must seed before build.

WHY THIS TEST EXISTS. `trainer.train()` calls `set_seed`, but it does so after its caller has
already constructed the network. So for the entire study `torch.manual_seed(7)` governed only
dropout and batch order, and every one of 945 deep-model fold-runs drew its weights from an
unseeded RNG. Nothing failed, because every assertion in this project checks a VALUE and none
checked whether the value was re-derivable.

The panel means were not threatened -- per-dataset training noise of 0.006-0.010 induces about
0.001 on a mean over 94 datasets -- but exact reproducibility was gone and the manuscript's
"identical seed" was false.

Two tests: the behaviour (same seed gives the same weights) and the ordering (the scripts seed
before they build), because the behaviour test passes even in the broken code if the caller
happens to seed first, and the ordering test is what actually regressed.
"""

import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

torch = pytest.importorskip("torch")

from rbp.models import registry  # noqa: E402
from rbp.train import trainer  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402


def _first_weight_sum(model_name, seed):
    trainer.set_seed(seed)
    h = registry.build(model_name, cfgmod.load(ROOT / "config" / "params.yaml"))
    return float(next(p for p in h.model.parameters() if p.requires_grad).sum())


def test_same_seed_gives_identical_initial_weights():
    a = _first_weight_sum("cnn", 7)
    b = _first_weight_sum("cnn", 7)
    assert a == b, "seeding before build must make initialisation deterministic"


def test_different_seeds_give_different_initial_weights():
    """Guards the reverse failure: a test that passes because nothing is random at all."""
    assert _first_weight_sum("cnn", 7) != _first_weight_sum("cnn", 8)


@pytest.mark.parametrize("script", ["scripts/cloud_train.py", "scripts/train.py"])
def test_entry_points_seed_before_they_build(script):
    """The ordering IS the bug. Read the source; a value check cannot see this."""
    src = (ROOT / script).read_text()
    seed_at = [m.start() for m in re.finditer(r"trainer\.set_seed\(", src)]
    build_at = [m.start() for m in re.finditer(r"registry\.build\(", src)]
    assert seed_at, f"{script} never calls trainer.set_seed before building a model"
    assert build_at, f"{script} does not build a model; test needs updating"
    assert min(seed_at) < min(build_at), (
        f"{script} calls registry.build before trainer.set_seed, so model weights are "
        f"initialised from an unseeded RNG")
