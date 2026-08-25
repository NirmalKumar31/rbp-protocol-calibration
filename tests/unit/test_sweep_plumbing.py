"""The parts of the GPU sweep that decide where work goes and where results land.

No torch and no network: these are the pure functions inside scripts/cloud_train.py, and
they are the ones whose failure is silent. A run prefix missing an axis overwrites another
run's results; a broken stride either drops manifest rows or runs them twice, and in both
cases every individual task still exits 0.
"""

import importlib.util
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))


def _load():
    """Import scripts/cloud_train.py by path: scripts/ is not a package."""
    spec = importlib.util.spec_from_file_location(
        "cloud_train", ROOT / "scripts" / "cloud_train.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ct = _load()

AXES = dict(arm="dinuc", cell="K562", protein="QKI", model="cnn", fold=0)


class TestRunPrefix:
    def test_contains_every_axis(self):
        p = ct.run_prefix(**AXES)
        for value in ("dinuc", "K562", "QKI", "cnn", "fold0"):
            assert value in p

    @pytest.mark.parametrize("axis,other", [
        ("arm", "gc"), ("cell", "HepG2"), ("protein", "PUM2"),
        ("model", "splicebert"), ("fold", 3)])
    def test_changing_any_axis_changes_the_prefix(self, axis, other):
        """The collision bug, generalised. Three separate defects in this project were
        two runs that had to differ writing to one path: the panel files keyed only by
        cell line, the default outdir, and the five folds of one model."""
        assert ct.run_prefix(**AXES) != ct.run_prefix(**{**AXES, axis: other})

    def test_prefixes_are_unique_across_a_realistic_grid(self):
        seen = {ct.run_prefix("dinuc", c, p, m, f)
                for c in ("K562", "HepG2")
                for p in ("QKI", "PUM2", "BCLAF1")
                for m in ("cnn", "rnabert", "splicebert")
                for f in range(5)}
        assert len(seen) == 2 * 3 * 3 * 5


class TestTaskIndex:
    """SHARD/NSHARDS striding. Five Batch jobs share one manifest."""

    @staticmethod
    def _indices(nshards, counts):
        """What the whole fleet would resolve, given each shard's task count."""
        out = []
        for shard in range(nshards):
            for local in range(counts[shard]):
                os.environ.update(SHARD=str(shard), NSHARDS=str(nshards),
                                  BATCH_TASK_INDEX=str(local))
                out.append(ct.task_index(sum(counts)))
        return out

    def teardown_method(self):
        for k in ("SHARD", "NSHARDS", "BATCH_TASK_INDEX"):
            os.environ.pop(k, None)

    def test_single_shard_is_the_identity(self):
        assert self._indices(1, [7]) == list(range(7))

    def test_five_shards_partition_the_manifest_exactly_once(self):
        total, n = 4725, 5
        counts = [(total - i + n - 1) // n for i in range(n)]
        assert sum(counts) == total
        got = self._indices(n, counts)
        assert sorted(got) == list(range(total))

    def test_an_uneven_manifest_still_partitions(self):
        total, n = 13, 5
        counts = [(total - i + n - 1) // n for i in range(n)]
        assert sorted(self._indices(n, counts)) == list(range(total))

    def test_shards_interleave_rather_than_slice(self):
        """The manifest is ordered most-expensive first. A contiguous slice would hand
        one region every large run and leave it going hours after the rest drained."""
        os.environ.update(SHARD="0", NSHARDS="5", BATCH_TASK_INDEX="0")
        first = ct.task_index(100)
        os.environ["BATCH_TASK_INDEX"] = "1"
        second = ct.task_index(100)
        assert (first, second) == (0, 5)

    def test_defaults_to_zero_outside_batch(self):
        assert ct.task_index(10) == 0
