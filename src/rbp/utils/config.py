"""Load config/params.yaml and resolve paths relative to the project root."""

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[3]
DEFAULT = ROOT / "config" / "params.yaml"


class Config(dict):
    """Dict with attribute access and project-root-relative path helpers."""

    def __getattr__(self, k):
        try:
            v = self[k]
        except KeyError as e:
            raise AttributeError(k) from e
        return Config(v) if isinstance(v, dict) else v

    def path(self, key, *parts):
        """Resolve a `paths:` entry, creating the directory."""
        p = ROOT / self["paths"][key]
        for part in parts:
            p = p / str(part)
        (p.parent if p.suffix else p).mkdir(parents=True, exist_ok=True)
        return p


def load(path=None):
    path = Path(path) if path else DEFAULT
    with open(path) as fh:
        return Config(yaml.safe_load(fh))


def proteins(cfg=None, panel=None):
    """Protein panel as a list of dicts. Defaults to config/proteins.tsv.

    `panel` selects a different TSV with the same columns, which is how the full
    139-protein sweep runs alongside the development panel without either
    overwriting the other.
    """
    f = Path(panel) if panel else ROOT / "config" / "proteins.tsv"
    if not f.is_absolute():
        f = ROOT / f
    if not f.exists():
        raise FileNotFoundError(f"{f} missing - run `make panel` first")
    lines = f.read_text().strip().splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, ln.split("\t"))) for ln in lines[1:]]


def protein_names(cfg=None):
    return [p["protein"] for p in proteins(cfg)]
