"""Where the cloud lives. One source of truth for project and bucket names.

WHY THIS FILE EXISTS. Eighteen files hardcoded the string "rbp-composition-2026". That is
invisible for as long as you only ever run in that project, and it is a total failure the
first time somebody tries to reproduce the work somewhere else -- which is the whole point
of a reproducible pipeline. A hardcoded project id is not a small tidiness problem, it is
the difference between "runs anywhere" and "runs on the author's account".

RESOLUTION ORDER, most specific first:

    1. the explicit argument, if a caller passes one
    2. the environment: GOOGLE_CLOUD_PROJECT / DERIVED_BUCKET / RAW_BUCKET
    3. config/params.yaml -> cloud:
    4. nothing. It raises.

There is deliberately NO fallback default. A default is how you end up writing results into
somebody else's bucket, or reading a stale one and never noticing the run did nothing. If
the environment is not configured, that is a bug in the environment and it should stop.

BUCKET NAMING. Buckets are globally unique across all of Google Cloud, so a reproducer
cannot have `rbp-composition-2026-derived`; someone already does. The convention is
`{project_id}-derived` and `{project_id}-raw`, which inherits uniqueness from the project id
and means one variable configures everything.
"""

import os
from functools import lru_cache

ENV_PROJECT = "GOOGLE_CLOUD_PROJECT"
ENV_DERIVED = "DERIVED_BUCKET"
ENV_RAW = "RAW_BUCKET"


@lru_cache(maxsize=1)
def _from_config():
    """The `cloud:` block of params.yaml, or {} if the file has none."""
    try:
        from . import config as cfgmod
        cfg = cfgmod.load()
        return dict(cfg["cloud"]) if "cloud" in cfg else {}
    except Exception:                       # config is optional for cloud resolution
        return {}


def _resolve(explicit, env_key, cfg_key, what):
    if explicit:
        return explicit
    if os.environ.get(env_key):
        return os.environ[env_key]
    v = _from_config().get(cfg_key)
    if v:
        return v
    raise RuntimeError(
        f"{what} is not configured. Set ${env_key}, or add cloud.{cfg_key} to "
        f"config/params.yaml. There is no default on purpose -- guessing a project or "
        f"bucket name is how a run silently reads or writes the wrong account.")


def project(explicit=None):
    return _resolve(explicit, ENV_PROJECT, "project_id", "GCP project")


def derived_bucket(explicit=None):
    """Derived artefacts: processed datasets, panels, manifests, runs, results."""
    if explicit:
        return explicit
    if os.environ.get(ENV_DERIVED):
        return os.environ[ENV_DERIVED]
    v = _from_config().get("derived_bucket")
    return v or f"{project()}-derived"


def raw_bucket(explicit=None):
    """Immutable inputs: genome, annotation, ClinVar, ENCODE peaks."""
    if explicit:
        return explicit
    if os.environ.get(ENV_RAW):
        return os.environ[ENV_RAW]
    v = _from_config().get("raw_bucket")
    return v or f"{project()}-raw"


def client(explicit_project=None):
    """A storage client bound to the resolved project."""
    from google.cloud import storage
    return storage.Client(project=project(explicit_project))


def bucket(name=None, explicit_project=None):
    return client(explicit_project).bucket(derived_bucket(name))


def describe():
    """One line for a log header, so every run records where it was pointed."""
    try:
        return f"project={project()} derived={derived_bucket()} raw={raw_bucket()}"
    except RuntimeError as e:
        return f"UNCONFIGURED: {e}"
