"""A directory that answers to the google-cloud-storage Bucket interface.

WHY THIS EXISTS. The GCP project's billing account was closed, so
`gs://rbp-repro-2026-derived` returns 403 on every object and the sweep can no longer read
its own inputs. Every input still needed is on local disk, so the missing piece is not the
data but the API: `cloud_train.py` talks to GCS through a small, well-defined surface --
eight blob methods and one list call -- and nothing else in it knows about buckets.

So rather than rewrite the sweep, implement that surface over a directory. `cloud_train.py`
changes in exactly one function body; its manifest, completion-marker, checkpoint and
aggregate logic are untouched, which is what lets the GC arm claim the same protocol as the
dinucleotide arm that has already run.

    store = LocalBucket("/store")           # writes and reads under /store
    store = LocalBucket("/tmp/out", ro="/vol")   # writes to /tmp/out, reads fall through

The read-through root is how a Modal task works: the datasets arrive on a Volume mounted
read-only, and the run's outputs go to container-local disk to be returned to the driver.
Nothing writes to shared storage during a sweep, so there is no commit to conflict over.

WRITES ARE ATOMIC, and that is load-bearing rather than tidiness. The sweep's resume rule is
"a run whose metrics.json exists is done". GCS gives that for free because an object appears
whole or not at all. On a filesystem a task killed mid-write leaves a truncated marker that
reads as complete, so a lost run would be silently skipped and its dataset would carry four
folds instead of five. Write to a temporary name in the same directory, then rename.
"""

import os
import shutil
from pathlib import Path


class NotFound(Exception):
    """Raised by Blob.delete on a missing object, mirroring google.api_core."""


class Blob:
    def __init__(self, root, name, ro=None):
        self.name = name
        self._w = Path(root) / name
        self._ro = Path(ro) / name if ro else None

    @property
    def path(self):
        """Where this blob actually is: the writable copy wins, else the read-only root."""
        if self._w.exists() or self._ro is None:
            return self._w
        return self._ro if self._ro.exists() else self._w

    def exists(self):
        return self.path.exists()

    def download_as_bytes(self):
        return self.path.read_bytes()

    def download_as_text(self):
        return self.path.read_text()

    def download_to_filename(self, dst):
        dst = Path(dst)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(self.path, dst)

    def upload_from_string(self, body, content_type=None):
        self._write(body.encode() if isinstance(body, str) else body)

    def upload_from_filename(self, src):
        self._write(Path(src).read_bytes())

    def delete(self):
        if not self._w.exists():
            raise NotFound(self.name)
        self._w.unlink()

    def _write(self, data):
        self._w.parent.mkdir(parents=True, exist_ok=True)
        tmp = self._w.with_name(self._w.name + f".tmp{os.getpid()}")
        tmp.write_bytes(data)
        os.replace(tmp, self._w)          # atomic within a filesystem


class _Client:
    """Only list_blobs is ever called on the client, and it ignores the bucket name."""

    def __init__(self, root, ro=None):
        self.root, self.ro = Path(root), Path(ro) if ro else None

    def list_blobs(self, _bucket=None, prefix=""):
        seen = set()
        for base in (self.root, self.ro):
            if base is None or not base.exists():
                continue
            start = base / prefix
            if not start.exists():
                continue
            for p in sorted(start.rglob("*")):
                if not p.is_file():
                    continue
                name = str(p.relative_to(base))
                if name in seen or name.endswith(".tmp") or ".tmp" in p.name:
                    continue
                seen.add(name)
                yield Blob(self.root, name, self.ro)


class LocalBucket:
    def __init__(self, root, ro=None):
        self.root = Path(root)
        self.ro = Path(ro) if ro else None
        self.name = str(self.root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.client = _Client(self.root, self.ro)

    def blob(self, name):
        return Blob(self.root, name, self.ro)


def uri(bucket, key=""):
    """Name an object for a log line without claiming gs:// for a directory."""
    base = bucket.name if isinstance(bucket, LocalBucket) else f"gs://{bucket.name}"
    return f"{base}/{key}" if key else base


def from_env():
    """A LocalBucket if STORE_DIR is set, else None so the caller falls back to GCS."""
    d = os.environ.get("STORE_DIR")
    return LocalBucket(d, os.environ.get("STORE_RO")) if d else None
