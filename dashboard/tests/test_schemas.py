"""The tables the figures read must exist, carry their columns, and never be written to.

These are cheap checks against an expensive failure: a renamed column upstream would otherwise
surface as a traceback inside a chart, in a browser, in front of whoever the dashboard was built
to impress.
"""

from __future__ import annotations

import ast
import csv
from pathlib import Path

import pytest
from rbp_dashboard import copy as dash_copy
from rbp_dashboard import data, figures, theme

PACKAGE = Path(data.__file__).parent


@pytest.mark.parametrize("name", sorted(data.REQUIRED))
def test_every_required_table_exists_and_is_not_empty(name):
    path = data.TABLES / name
    assert path.exists(), f"{name} is missing from results/tables/"
    with path.open() as handle:
        rows = list(csv.reader(handle))
    assert len(rows) >= 2, f"{name} has no data rows"


@pytest.mark.parametrize("name,columns", sorted(data.REQUIRED.items()))
def test_required_columns_are_present(name, columns):
    with (data.TABLES / name).open() as handle:
        header = next(csv.reader(handle))
    missing = [c for c in columns if c not in header]
    assert not missing, f"{name} lost {missing}"


def test_load_raises_schema_error_for_a_missing_table():
    with pytest.raises(data.SchemaError):
        data.load("this_table_does_not_exist.csv")


def test_check_refuses_an_ambiguous_or_absent_label():
    with pytest.raises(data.SchemaError):
        data.check("cross_fitting.csv", "a label that is not in the table")


def test_no_dashboard_module_can_write_to_the_repository():
    """Parsed, not trusted. `open(..., "w")`, `to_csv`, `unlink` and friends must not appear.

    The dashboard reads published evidence. A stray write here would corrupt the artefact the
    paper cites, so this is checked structurally rather than left to review.
    """
    banned_calls = {
        "to_csv", "to_parquet", "to_json", "to_excel", "to_feather", "to_pickle",
        "write", "writelines", "writerow", "writerows", "unlink", "rmtree", "remove",
        "rename", "replace", "mkdir", "touch", "write_text", "write_bytes", "savefig",
        "write_image", "write_html", "chmod",
    }
    offenders = []
    for path in sorted(PACKAGE.glob("*.py")) + [PACKAGE.parent / "app.py"]:
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                func = node.func
                name = getattr(func, "attr", None) or getattr(func, "id", None)
                if name in banned_calls:
                    offenders.append(f"{path.name}:{node.lineno} {name}()")
                if name == "open":
                    for arg in list(node.args[1:]) + [
                        k.value for k in node.keywords if k.arg == "mode"
                    ]:
                        if isinstance(arg, ast.Constant) and any(
                            m in str(arg.value) for m in ("w", "a", "x", "+")
                        ):
                            offenders.append(f"{path.name}:{node.lineno} open(mode={arg.value!r})")
    assert not offenders, "dashboard code contains write operations:\n  " + "\n  ".join(offenders)


def test_the_palette_is_a_fixed_order_of_three_distinct_colours():
    assert len(theme.CATEGORICAL) == 3
    assert len(set(theme.CATEGORICAL)) == 3
    assert set(theme.PROTOCOL_COLOUR) == {"gc", "dn", "neg2"}
    assert set(theme.MODEL_COLOUR) == {"kmer", "cnn", "splicebert"}
    # Not purple, which was an explicit requirement. Reject hues where blue and red both
    # dominate green, which is what reads as purple or violet.
    for colour in theme.CATEGORICAL:
        r, g, b = (int(colour[i:i + 2], 16) for i in (1, 3, 5))
        assert not (b > g + 40 and r > g + 40), f"{colour} reads as purple"


def test_every_figure_builds_from_the_committed_tables():
    """Each figure is constructed once. A figure that raises is a broken view."""
    contrast = data.table("three_arm_contrast.csv")
    per_dataset = data.table("three_arm_per_dataset.csv")
    models = data.table("three_arm_models.csv")
    cross = data.table("cross_fitting.csv")
    external = data.table("external_replication.csv")

    built = [
        figures.protocol_levels(contrast),
        figures.protocol_contribution(contrast),
        figures.contribution_distribution(per_dataset),
        figures.paired_slopes(per_dataset, "gc", "dn"),
        figures.apparent_versus_contribution(per_dataset),
        figures.model_spans(models),
        figures.model_by_protocol(models),
        figures.known_null(cross),
        figures.estimator_comparison(cross),
        figures.external_arms(data.table("external_replication_per_dataset.csv")),
        figures.external_levels(external, "baseline"),
        figures.external_levels(external, "contribution"),
        figures.provenance_mix(data.table("PROVENANCE.csv")),
    ]
    assert all(fig.data for fig in built), "a figure was built with no traces"


def test_no_figure_uses_a_second_y_axis():
    """One axis, always. Two scales in one frame is the most common chart defect there is."""
    contrast = data.table("three_arm_contrast.csv")
    for fig in (figures.protocol_levels(contrast), figures.protocol_contribution(contrast)):
        assert "yaxis2" not in fig.layout, "a figure grew a second y-axis"


def test_every_glossary_term_used_in_prose_is_defined():
    for key in ("span", "nested contribution", "cross-fitted estimator", "two-stage estimator",
                "held-out benchmark", "GC-matched", "dinucleotide-matched", "bias-aware",
                "estimator floor"):
        assert key in dash_copy.GLOSSARY, f"{key} is referenced but not defined"
        assert len(dash_copy.GLOSSARY[key]) > 40, f"{key} has a stub definition"
