"""Read-only access to the committed tables in `results/tables/`.

Every function here opens a file for reading and returns an in-memory frame. Nothing in this
module writes, moves or deletes anything, and `tests/test_schemas.py` asserts that no module in
this package contains a write call. Filtering happens on copies.

Two things this module refuses to do. It does not compute a headline estimate, and it does not
hold one as a literal. Each headline is pulled from the `check` row that the release already
asserts, so the dashboard cannot disagree with the paper without a test failing.
"""

from __future__ import annotations

import functools
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
TABLES = REPO_ROOT / "results" / "tables"

# Column sets checked before a view renders. Chosen as the columns the figures actually read,
# not the full header, so an added column upstream does not fail the dashboard.
REQUIRED = {
    "three_arm_contrast.csv": ["check", "value", "ci_low", "ci_high", "n"],
    "three_arm_per_dataset.csv": [
        "dataset", "protein", "cell",
        "comp_gc", "full_gc", "gain_gc", "n_gc",
        "comp_dn", "full_dn", "gain_dn", "n_dn",
        "comp_neg2", "full_neg2", "gain_neg2", "n_neg2",
    ],
    "three_arm_models.csv": ["check", "value", "ci_low", "ci_high", "n"],
    "three_arm_models_per_dataset.csv": [
        "dataset", "protein", "cell",
        "kmer_gain_gc", "cnn_gain_gc", "splicebert_gain_gc",
        "kmer_gain_dn", "cnn_gain_dn", "splicebert_gain_dn",
        "kmer_gain_neg2", "cnn_gain_neg2", "splicebert_gain_neg2",
    ],
    "cross_fitting.csv": ["check", "value", "ci_low", "ci_high", "n"],
    "estimator_floor.csv": ["check", "value", "ci_low", "ci_high", "n"],
    "external_replication.csv": ["check", "value", "ci_low", "ci_high", "n"],
    "external_replication_per_dataset.csv": [
        "dataset", "protein", "cell", "comp_n1", "gain_n1", "comp_n2", "gain_n2",
    ],
    "panel_summary.csv": ["dataset", "protein", "cell", "pairs", "in_both_arms"],
    "PROVENANCE.csv": ["table", "producing_script", "status", "sha256", "bytes"],
    "match_quality.csv": ["check", "value", "n"],
    "cost_of_matching.csv": ["dataset", "protein", "cell", "pairs", "cost",
                             "auroc_gc", "auroc_dn"],
    "negative_draws.csv": ["check", "value", "n"],
    "class_ratio.csv": ["check", "value", "n"],
    "negative_set_survey.csv": ["check", "value", "n"],
    "negative_set_survey_per_method.csv": ["method", "year", "kind",
                                           "composition_baseline", "url", "quote"],
    "recommendation_works.csv": ["check", "value", "ci_low", "ci_high"],
    "variant_ladder.csv": ["arm", "n_variants", "auroc"],
}


class SchemaError(RuntimeError):
    """A table is missing, empty, or lacks a column a figure needs."""


@functools.cache
def load(name: str) -> pd.DataFrame:
    path = TABLES / name
    if not path.exists():
        raise SchemaError(f"{name} is not in results/tables/")
    # `float_precision="round_trip"` is not a stylistic choice. pandas' default C parser is
    # faster and slightly lossy: it reads 0.06626403675115053 as 0.0662640367511505, a
    # difference of about 1e-17. That is far below anything displayed, but it means the
    # dashboard would be showing a re-parsed approximation of a published number rather than
    # the number. The round-trip parser reproduces the committed digits exactly, and
    # tests/test_headline_values.py compares against the stdlib float to keep it that way.
    frame = pd.read_csv(path, float_precision="round_trip")
    if frame.empty:
        raise SchemaError(f"{name} has a header but no rows")
    missing = [c for c in REQUIRED.get(name, []) if c not in frame.columns]
    if missing:
        raise SchemaError(f"{name} is missing {', '.join(missing)}")
    return frame


def table(name: str) -> pd.DataFrame:
    """A copy, so a caller that filters cannot mutate the cached frame."""
    return load(name).copy()


def check(name: str, label: str) -> pd.Series:
    """One asserted row, matched on the exact `check` string the release uses.

    Exact rather than substring: `"nested contribution, gc arm"` is a prefix of nothing here,
    but `"4-mer contribution as published, gc arm"` and `"...fully cross-fitted, gc arm"` share
    a long prefix, and a substring match would silently pick whichever came first.
    """
    frame = load(name)
    hit = frame[frame["check"] == label]
    if len(hit) != 1:
        raise SchemaError(f"{name}: expected exactly one row for {label!r}, found {len(hit)}")
    return hit.iloc[0]


def value(name: str, label: str) -> float:
    return float(check(name, label)["value"])


def interval(name: str, label: str) -> tuple[float, float | None, float | None]:
    row = check(name, label)
    low, high = row.get("ci_low"), row.get("ci_high")
    return (
        float(row["value"]),
        None if pd.isna(low) else float(low),
        None if pd.isna(high) else float(high),
    )


def ci_text(name: str, label: str, digits: int = 3) -> str:
    val, low, high = interval(name, label)
    if low is None or high is None:
        return f"{val:.{digits}f}"
    return f"{val:.{digits}f} (95% CI {low:.{digits}f} to {high:.{digits}f})"


def available() -> list[str]:
    """Which required tables are present. Used by the app to fail with a list, not a stack."""
    return [n for n in REQUIRED if (TABLES / n).exists()]


# --------------------------------------------------------------------------- interactive filters

def facets(frame: pd.DataFrame) -> dict:
    """The filter options a per-dataset table can offer, read off the table itself.

    Read rather than hardcoded, so a regenerated panel with a third cell line offers it without
    an edit here.
    """
    return {
        "cells": sorted(frame["cell"].dropna().unique().tolist()),
        "proteins": sorted(frame["protein"].dropna().unique().tolist()),
    }


def apply_filters(
    frame: pd.DataFrame,
    cells: list[str] | None = None,
    proteins: list[str] | None = None,
    size_range: tuple[float, float] | None = None,
    size_column: str | None = None,
) -> pd.DataFrame:
    """Filter a per-dataset table. Returns a copy; the cached frame is never touched."""
    out = frame.copy()
    if cells:
        out = out[out["cell"].isin(cells)]
    if proteins:
        out = out[out["protein"].isin(proteins)]
    if size_range and size_column and size_column in out.columns:
        low, high = size_range
        out = out[(out[size_column] >= low) & (out[size_column] <= high)]
    return out


def descriptive(frame: pd.DataFrame, column: str) -> dict | None:
    """Summary statistics of a FILTERED subset. Descriptive only, never an estimand.

    This is the one place the dashboard produces a number that is not in a committed table, and
    the distinction matters enough to state twice. A published contribution is a panel mean over
    all 94 datasets with a protein-clustered bootstrap interval. What this returns is the mean
    of whatever the user filtered to, with no interval, no clustering and no claim. The interface
    marks every one of these with a LIVE badge and every published value with a PUBLISHED badge.

    `interval` is deliberately absent from the return value so a caller cannot draw an error bar
    on it by reaching for a key that happens to exist.
    """
    series = frame[column].dropna()
    if series.empty:
        return None
    return {
        "n": int(series.size),
        "mean": float(series.mean()),
        "median": float(series.median()),
        "min": float(series.min()),
        "max": float(series.max()),
    }
