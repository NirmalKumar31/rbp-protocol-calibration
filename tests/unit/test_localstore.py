"""The local store must behave like a GCS bucket, including where it is easy not to.

`cloud_train.py` was written against google-cloud-storage and its resume rule leans on two
properties of that API which a naive directory does NOT have:

  * an object appears whole or not at all, so a completion marker is never half-written;
  * deleting something absent raises NotFound rather than succeeding silently.

Both matter. The sweep skips any run whose metrics.json exists, so a truncated marker means
a dataset silently carries four folds instead of five -- and every downstream AUROC is then
computed on 80% of the rows with nothing anywhere saying so.
"""

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from rbp.utils.localstore import LocalBucket, NotFound, uri  # noqa: E402


def test_round_trip_text_and_bytes(tmp_path):
    b = LocalBucket(tmp_path)
    b.blob("a/b.txt").upload_from_string("hello")
    b.blob("a/c.bin").upload_from_string(b"\x00\x01")
    assert b.blob("a/b.txt").download_as_text() == "hello"
    assert b.blob("a/c.bin").download_as_bytes() == b"\x00\x01"


def test_missing_blob_does_not_exist(tmp_path):
    assert not LocalBucket(tmp_path).blob("nope").exists()


def test_upload_and_download_filename(tmp_path):
    src = tmp_path / "src.txt"
    src.write_text("payload")
    b = LocalBucket(tmp_path / "store")
    b.blob("x/y.txt").upload_from_filename(src)
    out = tmp_path / "deep/out.txt"
    b.blob("x/y.txt").download_to_filename(out)
    assert out.read_text() == "payload"


def test_delete_missing_raises_not_found(tmp_path):
    b = LocalBucket(tmp_path)
    b.blob("gone").upload_from_string("x")
    b.blob("gone").delete()
    with pytest.raises(NotFound):
        b.blob("gone").delete()


def test_writes_are_atomic_and_leave_no_temp_files(tmp_path):
    b = LocalBucket(tmp_path)
    for i in range(5):
        b.blob(f"runs/f{i}/metrics.json").upload_from_string(json.dumps({"i": i}))
    leftovers = [p for p in Path(tmp_path).rglob("*") if ".tmp" in p.name]
    assert leftovers == []


def test_list_blobs_respects_prefix(tmp_path):
    b = LocalBucket(tmp_path)
    b.blob("runs/gc/a/metrics.json").upload_from_string("{}")
    b.blob("runs/gc/b/metrics.json").upload_from_string("{}")
    b.blob("other/thing").upload_from_string("x")
    names = sorted(x.name for x in b.client.list_blobs("bucket", prefix="runs/gc/"))
    assert names == ["runs/gc/a/metrics.json", "runs/gc/b/metrics.json"]


def test_read_through_falls_back_to_readonly_root(tmp_path):
    ro, rw = tmp_path / "vol", tmp_path / "out"
    LocalBucket(ro).blob("processed/gc/K562/X/dataset.tsv").upload_from_string("data")
    b = LocalBucket(rw, ro=ro)
    assert b.blob("processed/gc/K562/X/dataset.tsv").download_as_text() == "data"


def test_writable_root_shadows_readonly(tmp_path):
    ro, rw = tmp_path / "vol", tmp_path / "out"
    LocalBucket(ro).blob("k").upload_from_string("old")
    b = LocalBucket(rw, ro=ro)
    b.blob("k").upload_from_string("new")
    assert b.blob("k").download_as_text() == "new"
    assert (ro / "k").read_text() == "old"      # the read-only root is never touched


def test_list_blobs_spans_both_roots_without_duplicates(tmp_path):
    ro, rw = tmp_path / "vol", tmp_path / "out"
    LocalBucket(ro).blob("runs/a").upload_from_string("x")
    b = LocalBucket(rw, ro=ro)
    b.blob("runs/a").upload_from_string("y")
    b.blob("runs/b").upload_from_string("z")
    names = sorted(x.name for x in b.client.list_blobs(prefix="runs"))
    assert names == ["runs/a", "runs/b"]


def test_uri_never_claims_gs_for_a_directory(tmp_path):
    assert uri(LocalBucket(tmp_path), "k").startswith(str(tmp_path))
    assert "gs://" not in uri(LocalBucket(tmp_path))
