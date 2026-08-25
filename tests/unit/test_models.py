"""Model encoding, architecture, config wiring, and the checkpoint/resume contract.

The language models are not built here: instantiating them downloads hundreds of MB.
Their integration is covered by the smoke test instead. What is tested here is
everything that can go wrong without a network.
"""

import json

import numpy as np
import pytest
import torch

from rbp.models import cnn as cnn_mod
from rbp.models import registry
from rbp.train import data as tdata
from rbp.train import trainer


class TestOneHot:
    def test_shape_and_values(self):
        x = cnn_mod.one_hot("ACGU")
        assert x.shape == (4, 4)
        assert torch.equal(x, torch.eye(4))

    def test_one_hot_per_position(self):
        x = cnn_mod.one_hot("AAGGUUCC")
        assert torch.equal(x.sum(dim=0), torch.ones(8))

    def test_unknown_base_is_all_zero(self):
        x = cnn_mod.one_hot("ANA")
        assert x[:, 1].sum() == 0
        assert x[:, 0].sum() == 1

    def test_thymine_is_not_recognised(self):
        # windows are stored as RNA; a stray T must not silently encode as U
        assert cnn_mod.one_hot("T").sum() == 0

    def test_batch_shape(self):
        x = cnn_mod.one_hot_batch(["ACGU", "GGGG"])
        assert x.shape == (2, 4, 4)


class TestCNN:
    def test_forward_shape(self):
        m = cnn_mod.DeepBindCNN()
        out = m(cnn_mod.one_hot_batch(["ACGU" * 25 + "A"] * 3))
        assert out.shape == (3,)

    def test_param_count_is_small(self):
        # the CNN is meant to be the tiny honest floor, not a big model
        assert 5_000 < cnn_mod.DeepBindCNN().n_params < 12_000

    def test_global_pool_makes_position_nearly_irrelevant(self):
        """Moving a motif around the window should barely change the score.

        The EDA showed the signal sits ~15nt off centre by a protein-specific amount,
        so the model must not care much where a pattern falls. AdaptiveMaxPool1d(1)
        provides that.

        Invariance is approximate, not exact: MaxPool1d has stride 4, and the ReLU
        responses from the surrounding filler interact with bin boundaries, so shifts
        that are not clean multiples of the stride move the score slightly. The
        guarantee we rely on is that the variation is small, which this pins down.
        """
        m = cnn_mod.DeepBindCNN().eval()
        motif, filler = "GCAUG", "A"

        def score(pos, L=101):
            s = filler * pos + motif + filler * (L - pos - len(motif))
            with torch.no_grad():
                return float(m(cnn_mod.one_hot_batch([s])))

        interior = [score(p) for p in (20, 30, 40, 50, 60, 70)]
        assert np.ptp(interior) < 0.05          # untrained scale is ~0.1, so this is small

        # the structural guarantee: features are pooled to a single position, so the
        # classifier cannot see where in the window a pattern occurred
        with torch.no_grad():
            feats = m.features(cnn_mod.one_hot_batch([filler * 20 + motif + filler * 76]))
        assert feats.shape[-1] == 1

    def test_variable_length_is_accepted(self):
        m = cnn_mod.DeepBindCNN().eval()
        with torch.no_grad():
            assert m(cnn_mod.one_hot_batch(["ACGU" * 10])).shape == (1,)
            assert m(cnn_mod.one_hot_batch(["ACGU" * 40])).shape == (1,)


class TestRegistryConfig:
    """Model definitions must come from config, not be duplicated in code."""

    def test_names_come_from_config(self):
        names = registry.names()
        assert "cnn" in names
        assert len(names) >= 4

    def test_lora_threshold_decides_mode(self):
        cfg = {"models": {"big": {"kind": "lm", "params_m": 99.5},
                          "small": {"kind": "lm", "params_m": 0.5},
                          "tiny": {"kind": "cnn"}},
               "lora_threshold_m": 50}
        assert registry.mode("big", cfg) == "lora"
        assert registry.mode("small", cfg) == "full"
        assert registry.mode("tiny", cfg) == "scratch"

    def test_explicit_mode_overrides_threshold(self):
        cfg = {"models": {"m": {"kind": "lm", "params_m": 99.5, "mode": "full"}},
               "lora_threshold_m": 50}
        assert registry.mode("m", cfg) == "full"

    def test_unknown_model_raises(self):
        with pytest.raises(ValueError, match="unknown model"):
            registry.spec("not_a_model")

    def test_cnn_builds_from_config(self):
        h = registry.build("cnn")
        assert h.kind == "cnn" and h.mode == "scratch"
        assert h.sizes()["params_trainable"] == h.sizes()["params_total"]


class TestDataLoaders:
    def test_collate_returns_strings_and_float_labels(self):
        seqs, y = tdata.collate([("ACGU", 1.0), ("GGGG", 0.0)])
        assert seqs == ["ACGU", "GGGG"]
        assert y.dtype == torch.float32

    def test_class_balance_on_a_subset(self):
        from torch.utils.data import DataLoader, Subset
        ds = tdata.WindowDataset(["A"] * 10, [1] * 5 + [0] * 5)
        dl = DataLoader(Subset(ds, [0, 1, 2]), batch_size=2, collate_fn=tdata.collate)
        assert tdata.class_balance(dl) == {"n": 3, "positives": 3, "frac_positive": 1.0}


class TestCheckpointResume:
    """The contract that makes a preempted sweep cheap to restart."""

    def _loaders(self, n=64):
        rng = np.random.default_rng(0)
        # a learnable toy task: positives contain the motif
        pos = ["".join(rng.choice(list("ACGU"), 40)) + "GCAUG" for _ in range(n // 2)]
        neg = ["".join(rng.choice(list("ACGU"), 45)) for _ in range(n // 2)]
        seqs = pos + neg
        y = [1] * (n // 2) + [0] * (n // 2)
        ds = tdata.WindowDataset(seqs, y)
        from torch.utils.data import DataLoader
        dl = DataLoader(ds, batch_size=8, shuffle=True, collate_fn=tdata.collate)
        ev = DataLoader(ds, batch_size=8, collate_fn=tdata.collate)
        return {"train": dl, "val": ev, "test": ev}

    def test_checkpoint_roundtrip(self, tmp_path):
        h = registry.build("cnn")
        opt = torch.optim.AdamW(h.param_groups())
        ck = trainer.Checkpoint(tmp_path / "c.pt")
        assert not ck.exists()
        ck.save(h, opt, epoch=3, best=0.8, best_epoch=3,
                history=[{"epoch": 3}], elapsed=12.5)
        assert ck.exists()
        blob = ck.load()
        assert blob["epoch"] == 3 and blob["best"] == 0.8 and blob["elapsed"] == 12.5
        # best weights live in best.pt, not in here: they change only on improvement
        # while this blob is rewritten every epoch.
        assert "best_state" not in blob
        assert h.load(blob["model"]) == []          # no missing keys

    def test_finished_run_is_skipped(self, tmp_path):
        (tmp_path / "metrics.json").write_text(json.dumps({"test_auroc": 0.9}))
        h = registry.build("cnn")
        out = trainer.train(h, self._loaders(), tmp_path, epochs=99, log=lambda *_: None)
        assert out["test_auroc"] == 0.9            # returned without training

    def test_run_writes_all_artifacts(self, tmp_path):
        h = registry.build("cnn")
        m = trainer.train(h, self._loaders(), tmp_path, epochs=2, log=lambda *_: None)
        assert (tmp_path / "metrics.json").exists()
        assert (tmp_path / "best.pt").exists()
        assert (tmp_path / "test_predictions.npz").exists()
        # the checkpoint is removed once the run completes, so nothing stale is carried
        assert not (tmp_path / "checkpoint.pt").exists()
        assert m["epochs_run"] == 2 and len(m["history"]) == 2

    def test_resume_continues_from_the_saved_epoch(self, tmp_path):
        loaders = self._loaders()
        h1 = registry.build("cnn")
        opt = torch.optim.AdamW(h1.param_groups())
        trainer.Checkpoint(tmp_path / "checkpoint.pt").save(
            h1, opt, epoch=2, best=0.55, best_epoch=2,
            history=[{"epoch": 1, "val_auroc": 0.5, "train_loss": 0.7, "seconds": 1},
                     {"epoch": 2, "val_auroc": 0.55, "train_loss": 0.7, "seconds": 2}],
            elapsed=2.0)
        h2 = registry.build("cnn")
        m = trainer.train(h2, loaders, tmp_path, epochs=4, log=lambda *_: None)
        epochs = [r["epoch"] for r in m["history"]]
        assert epochs == [1, 2, 3, 4]              # picked up at 3, did not restart


    def test_resume_when_best_is_unrounded(self, tmp_path):
        """The bug that broke every preempted run.

        `best` is kept at full precision while history rounds val_auroc to 5 dp, so
        recovering best_epoch by matching the two found nothing and max() raised on an
        empty sequence. best_epoch is now persisted, so an unrounded best is harmless.
        """
        h = registry.build("cnn")
        opt = torch.optim.AdamW(h.param_groups())
        trainer.Checkpoint(tmp_path / "checkpoint.pt").save(
            h, opt, epoch=2, best=0.7973511219, best_epoch=2,
            history=[{"epoch": 1, "val_auroc": 0.79123, "train_loss": 0.7, "seconds": 1},
                     {"epoch": 2, "val_auroc": 0.79735, "train_loss": 0.7, "seconds": 2}],
            elapsed=2.0)
        m = trainer.train(registry.build("cnn"), self._loaders(), tmp_path,
                          epochs=3, log=lambda *_: None)
        assert [r["epoch"] for r in m["history"]] == [1, 2, 3]


class TestSweepSafety:
    """Everything the sweep needs that a laptop run does not.

    A spot preemption destroys the machine, so the only state that survives is what got
    mirrored off it. These pin the seam that makes that possible, and the two guards that
    stop a run producing a plausible number from a broken setup.
    """

    def _loaders(self, n=64):
        return TestCheckpointResume._loaders(self, n)

    def _one_class(self, split):
        dl = self._loaders()
        ds = tdata.WindowDataset(["ACGU" * 11] * 8, [1] * 8)
        from torch.utils.data import DataLoader
        dl[split] = DataLoader(ds, batch_size=4, collate_fn=tdata.collate)
        return dl

    # --- the on_epoch seam ------------------------------------------------------------
    def test_on_epoch_fires_once_per_epoch_after_the_checkpoint_exists(self, tmp_path):
        seen = []

        def hook(epoch, outdir):
            # the mirror must find a complete checkpoint, not one being written
            assert (outdir / trainer.CHECKPOINT).exists()
            seen.append(epoch)

        trainer.train(registry.build("cnn"), self._loaders(), tmp_path, epochs=3,
                      log=lambda *_: None, on_epoch=hook)
        assert seen == [1, 2, 3]

    def test_best_weights_are_on_disk_before_the_hook_runs(self, tmp_path):
        """best.pt is uploaded first. If it were written after the checkpoint, a kill
        between the two would leave a checkpoint naming weights that do not exist."""
        def hook(epoch, outdir):
            assert (outdir / trainer.BEST).exists()

        trainer.train(registry.build("cnn"), self._loaders(), tmp_path, epochs=2,
                      log=lambda *_: None, on_epoch=hook)

    def test_best_state_is_not_carried_in_the_checkpoint(self, tmp_path):
        """It changes only on improvement while the checkpoint is rewritten every epoch.
        For a full fine-tune it is a third of the payload."""
        holder = {}

        def hook(epoch, outdir):
            holder["blob"] = torch.load(outdir / trainer.CHECKPOINT, weights_only=False)

        trainer.train(registry.build("cnn"), self._loaders(), tmp_path, epochs=1,
                      log=lambda *_: None, on_epoch=hook)
        assert "best_state" not in holder["blob"]
        assert "elapsed" in holder["blob"]

    # --- timing across a resume -------------------------------------------------------
    def test_seconds_accumulates_across_a_resume(self, tmp_path):
        h = registry.build("cnn")
        opt = torch.optim.AdamW(h.param_groups())
        trainer.Checkpoint(tmp_path / trainer.CHECKPOINT).save(
            h, opt, epoch=1, best=0.5, best_epoch=1,
            history=[{"epoch": 1, "val_auroc": 0.5, "train_loss": 0.7, "seconds": 900.0}],
            elapsed=900.0)
        m = trainer.train(registry.build("cnn"), self._loaders(), tmp_path, epochs=2,
                          log=lambda *_: None)
        # without the fix this reported only the seconds since the resume, so a preempted
        # run under-reported its own cost by however long the first attempt ran
        assert m["seconds"] >= 900.0

    # --- degenerate splits ------------------------------------------------------------
    @pytest.mark.parametrize("split", ["train", "val", "test"])
    def test_a_one_class_split_is_refused_before_training(self, split, tmp_path):
        """nan AUROC never beats `best`, so best_state stays None, early stopping fires
        at `patience`, and the run reports val_auroc -1.0 with a test score from whatever
        the last epoch left behind. All of it silent."""
        with pytest.raises(ValueError, match="one class"):
            trainer.train(registry.build("cnn"), self._one_class(split), tmp_path,
                          epochs=2, log=lambda *_: None)

    def test_nothing_is_written_when_the_splits_are_refused(self, tmp_path):
        with pytest.raises(ValueError):
            trainer.train(registry.build("cnn"), self._one_class("val"), tmp_path,
                          epochs=2, log=lambda *_: None)
        assert not (tmp_path / "metrics.json").exists()

    # --- row ids ----------------------------------------------------------------------
    def test_ids_are_written_alongside_the_predictions(self, tmp_path):
        dl = self._loaders()
        ids = [f"r{i}" for i in range(len(dl["test"].dataset))]
        labels = trainer.labels_of(dl["test"])
        trainer.train(registry.build("cnn"), dl, tmp_path, epochs=1,
                      log=lambda *_: None, test_ids=(ids, labels))
        z = np.load(tmp_path / "test_predictions.npz", allow_pickle=True)
        assert list(z["id"]) == ids
        assert len(z["prob"]) == len(ids)

    def test_misaligned_ids_are_rejected(self, tmp_path):
        """A score vector attached to the wrong rows is the one error that would survive
        every downstream check: pooling, DeLong and the comparison against the
        composition arm would all run happily on scrambled data."""
        dl = self._loaders()
        labels = trainer.labels_of(dl["test"])
        ids = [f"r{i}" for i in range(len(labels))]
        with pytest.raises(ValueError, match="do not line up"):
            trainer.train(registry.build("cnn"), dl, tmp_path, epochs=1,
                          log=lambda *_: None, test_ids=(ids, labels[::-1]))

    def test_wrong_number_of_ids_is_rejected(self, tmp_path):
        dl = self._loaders()
        labels = trainer.labels_of(dl["test"])
        with pytest.raises(ValueError, match="do not line up"):
            trainer.train(registry.build("cnn"), dl, tmp_path, epochs=1,
                          log=lambda *_: None,
                          test_ids=(["r0", "r1"], labels[:2]))

    # --- ordering of the final writes --------------------------------------------------
    def test_metrics_json_is_written_after_the_predictions(self, tmp_path):
        """metrics.json is the completion marker locally and in GCS. It used to be
        written first, so a kill between the two left a run that looked finished and had
        no predictions."""
        trainer.train(registry.build("cnn"), self._loaders(), tmp_path, epochs=1,
                      log=lambda *_: None)
        assert ((tmp_path / "metrics.json").stat().st_mtime_ns
                >= (tmp_path / "test_predictions.npz").stat().st_mtime_ns)

    def test_skip_message_does_not_crash_on_a_shallow_outdir(self, tmp_path):
        """The log line used to index outdir.parents[2], which raises on a short path."""
        d = tmp_path / "run"
        d.mkdir()
        (d / "metrics.json").write_text(json.dumps({"test_auroc": 0.9}))
        out = trainer.train(registry.build("cnn"), self._loaders(), d, epochs=9,
                            log=lambda *_: None)
        assert out["test_auroc"] == 0.9
