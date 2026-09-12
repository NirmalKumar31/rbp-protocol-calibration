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


def test_the_decorative_nucleotide_palette_never_encodes_data():
    """NUCLEOTIDE fails CVD separation against itself and is a header motif, not a series.

    Green against cyan is dE 4.0 under tritanopia and green against rose is dE 4.6 under
    deuteranopia. That is fine for four stripes in a rule where nothing depends on telling them
    apart, and would be a defect in a legend. This asserts no figure reaches for it.
    """
    source = (PACKAGE / "figures.py").read_text()
    assert "NUCLEOTIDE" not in source, (
        "figures.py references theme.NUCLEOTIDE, which is decorative and not colour-vision safe "
        "as a categorical set")
    assert set(theme.NUCLEOTIDE) == {"A", "C", "G", "U"}


def test_the_dark_palette_is_distinct_and_readable_on_the_dark_surface():
    """The theme moved from light to dark, so the colours moved too. Check they actually did.

    An inverted light palette is the usual mistake. These are different hues chosen against
    SURFACE, and the surface must stay dark or the validation no longer describes what is shown.
    """
    r, g, b = (int(theme.SURFACE[i:i + 2], 16) for i in (1, 3, 5))
    assert r + g + b < 120, "SURFACE is no longer dark; the palette was validated against dark"
    for colour in theme.CATEGORICAL:
        cr, cg, cb = (int(colour[i:i + 2], 16) for i in (1, 3, 5))
        assert cr + cg + cb > 300, f"{colour} is too dark to read on {theme.SURFACE}"


def test_descriptive_statistics_carry_no_interval():
    """A filtered mean has no clustered interval. The return value must not offer one.

    If a `ci_low` ever appeared here, some figure would eventually draw an error bar on a
    filtered subset and present it as an estimate. Absence is the guard.
    """
    frame = data.table("three_arm_per_dataset.csv")
    stat = data.descriptive(data.apply_filters(frame, cells=["K562"]), "gain_dn")
    assert stat is not None
    assert set(stat) == {"n", "mean", "median", "min", "max"}
    for banned in ("ci_low", "ci_high", "interval", "ci", "se", "std"):
        assert banned not in stat


def test_filtering_never_mutates_the_cached_frame():
    before = len(data.load("three_arm_per_dataset.csv"))
    data.apply_filters(data.table("three_arm_per_dataset.csv"), cells=["K562"])
    assert len(data.load("three_arm_per_dataset.csv")) == before


def test_the_animated_figure_has_one_frame_per_model_class():
    fig = figures.animated_model_walk(data.table("three_arm_models.csv"))
    assert len(fig.frames) == 3
    assert {f.name for f in fig.frames} == {theme.MODEL_SHORT[m]
                                            for m in ("kmer", "cnn", "splicebert")}


def test_the_stylesheet_is_well_formed_and_its_font_url_is_on_one_line():
    """A CSS url() cannot contain a newline, and wrapping this one silently killed the fonts.

    It happened: the import was split across two source lines to satisfy the 100-character rule,
    which is valid Python, invalid CSS, and produces no error anywhere. The page simply fell back
    to a system font. Nothing else in the suite would have noticed.
    """
    import re

    css = theme.CSS
    assert css.lstrip().startswith("<style>") and css.rstrip().endswith("</style>")
    assert css.count("@import") == 1, "more than one font import, or none"
    assert "%(fonts)s" not in css, "the font placeholder was never substituted"
    match = re.search(r"@import url\('([^']+)'\)", css)
    assert match, "the font import is not a single well-formed url()"
    url = match.group(1)
    assert "\n" not in url and " " not in url, "the font URL is wrapped and will not load"
    assert "family=Inter" in url and "JetBrains+Mono" in url
    # The keyframe percentages must survive whatever assembly the module does to build this.
    assert "0%,100%" in css, "@keyframes percentages were eaten by string formatting"


def test_every_view_has_a_reading_guide_that_states_its_filter_behaviour():
    """The first build assumed the reader already knew the result. This is the correction.

    It also records whether the sidebar filter does anything on that view. Four of the seven show
    only panel-level published estimates, which correctly do not move when a protein is
    deselected, and saying nothing about that reads as a broken filter.
    """
    import importlib

    app = importlib.import_module("app")
    guides = set(dash_copy.GUIDE)
    views = set(app.VIEWS)
    assert views - guides == set(), f"views with no reading guide: {views - guides}"
    for name, guide in dash_copy.GUIDE.items():
        assert guide["lead"], f"{name} has no lead sentence"
        assert len(guide["checks"]) >= 3, f"{name} lists fewer than three things to check"
        assert guide["filters"] in dash_copy.FILTER_STATE, (
            f"{name} declares an unknown filter state {guide['filters']!r}")

    # Only views where the sidebar filter actually does something say so. Seven of ten were
    # printing "filters do not apply", which is a notice about a thing that is not happening.
    speaking = {k for k, (label, _) in dash_copy.FILTER_STATE.items() if label}
    assert speaking == {"responds", "partial"}, (
        f"only responsive views should carry a chip; {speaking} do")
    for state in speaking:
        assert dash_copy.FILTER_STATE[state][0] == "filters apply here", (
            "both speaking states should read the same, so the chip means one thing")
    assert dash_copy.FILTER_STATE["static"] == (None, None), (
        "a view the filter does not touch should say nothing at all")


def test_the_palette_is_muted_rather_than_saturated():
    """Separation is carried by lightness here, not chroma, and that is deliberate.

    Four desaturated sets were measured before this one and all four failed CVD separation:
    muting a hue removes chroma, and chroma is what the separation depends on. Spreading three
    colours across a wide lightness range keeps them calm and still measurably distinct. This
    checks the lightness spread is actually wide, so a future edit cannot quietly replace it
    with three colours of equal lightness that only pass on hue.
    """
    lightness = []
    for colour in theme.CATEGORICAL:
        r, g, b = (int(colour[i:i + 2], 16) / 255 for i in (1, 3, 5))
        lightness.append(0.2126 * r + 0.7152 * g + 0.0722 * b)
    assert max(lightness) - min(lightness) > 0.25, (
        f"the three colours are too close in lightness ({lightness}); this palette relies on "
        "lightness for separation because it is deliberately low in chroma")
    for colour in theme.CATEGORICAL:
        r, g, b = (int(colour[i:i + 2], 16) for i in (1, 3, 5))
        assert max(r, g, b) - min(r, g, b) < 130, f"{colour} is more saturated than this theme"


@pytest.mark.parametrize("width", [94, 49, 5, 2, 1, 0])
def test_every_view_renders_at_every_filter_width(width, monkeypatch):
    """Narrow filters are where widgets break, and the first harness never tested one.

    `st.slider("How many", 5, min(40, len(subset)))` raises when the subset has five or fewer
    rows, because min_value then meets or exceeds max_value and Streamlit refuses rather than
    clamping. It crashed in the browser at one dataset while the suite was green, because the
    only filters exercised left 49 and 94 rows. Zero rows is included: an empty filter is a
    thing a user does by accident and it must produce a message, not a traceback.
    """
    import importlib

    app = importlib.import_module("app")
    frame = data.table("three_arm_per_dataset.csv")
    proteins = sorted(frame["protein"].unique())

    if width == 94:
        flt = {"cells": None, "proteins": None, "size_range": None, "size_column": None}
    elif width == 49:
        flt = {"cells": ["K562"], "proteins": [], "size_range": None, "size_column": None}
    elif width == 0:
        flt = {"cells": ["no such cell line"], "proteins": [], "size_range": None,
               "size_column": None}
    else:
        flt = {"cells": None, "proteins": proteins[:width], "size_range": None,
               "size_column": None}

    for name, view in app.VIEWS.items():
        try:
            view(flt)
        except Exception as exc:  # noqa: BLE001
            pytest.fail(f"{name} raised at a filter width of {width}: "
                        f"{type(exc).__name__}: {exc}")


def test_every_chart_has_a_plain_language_caption_and_names_its_source():
    """Charts are drawn through `chart(fig, key, *sources)`. Both extras are required by shape.

    Parsed from the source rather than checked at runtime, so a caption missing from a branch
    that only renders under some filter still fails here. The audience for this dashboard
    includes people who do not work on RNA, and a chart with no sentence under it is not
    readable by them.
    """
    import ast

    app_src = (PACKAGE.parent / "app.py").read_text()
    tree = ast.parse(app_src)

    keys, figures_drawn = set(), set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        name = getattr(node.func, "id", None) or getattr(node.func, "attr", None)
        if name == "chart" and node.args:
            fig_arg = node.args[0]
            if isinstance(fig_arg, ast.Call):
                attr = getattr(fig_arg.func, "attr", None)
                if attr:
                    figures_drawn.add(attr)
            if len(node.args) > 1 and isinstance(node.args[1], ast.Constant):
                keys.add(node.args[1].value)

    assert keys, "no chart() calls found; the helper was bypassed"
    missing = sorted(k for k in keys if k not in dash_copy.PLAIN)
    assert not missing, f"charts drawn with no plain-language caption: {missing}"

    # The caption key is the figure function name. If they drift, a chart gets someone else's
    # sentence, which is worse than having none.
    mismatched = sorted(figures_drawn - set(dash_copy.PLAIN))
    assert not mismatched, f"figure functions with no caption entry: {mismatched}"

    for key, text in dash_copy.PLAIN.items():
        assert len(text) > 80, f"{key} has a stub caption"
        assert text.count("<b>") == text.count("</b>"), f"{key} has unbalanced markup"


def test_the_background_is_generated_and_decorative():
    """The helix field is computed from parameters and carries no data.

    Asserted because a background that looked like a plot would be a scientific claim. It is two
    sine waves and some rungs; nothing in it is measured, and it is not to scale.
    """
    svg = theme._helix_svg(800, 200, turns=4, amplitude=30, rungs=20)
    assert svg.startswith("<svg") and svg.endswith("</svg>")
    assert svg.count("<polyline") == 2, "a double helix has two strands"
    assert svg.count("<line") == 20, "rung count should follow the parameter"

    html = theme.background()
    assert 'aria-hidden="true"' in html, "decoration must be hidden from assistive technology"
    assert "prefers-reduced-motion" in theme.CSS, "the drift must be suppressible"
    # Layer count is asserted by test_the_background_is_deterministic_and_layered, which owns
    # it; duplicating the number here is how it came to say two after the field grew to four.


def test_the_schematic_diagrams_are_labelled_as_schematics_and_are_accessible():
    """Four method diagrams are drawn as page content. None is plotted from data.

    A diagram that looked like a plot would be a scientific claim, so each carries a visible
    SCHEMATIC label and a <title> that becomes its accessible name. They are content rather than
    decoration, so unlike the background they must NOT be aria-hidden.
    """
    from rbp_dashboard import graphics

    for name in ("hero", "negative_sets", "cross_fitting", "composition"):
        svg = getattr(graphics, name)()
        assert svg.startswith("<svg") and svg.rstrip().endswith("</svg>"), f"{name} is malformed"
        assert "<title>" in svg, f"{name} has no accessible name"
        assert "SCHEMATIC" in svg, f"{name} is not labelled as a schematic"
        assert 'aria-hidden' not in svg, f"{name} is content and must not be hidden"
        assert 'role="img"' in svg, f"{name} needs an image role"
        assert "<animate" in svg, f"{name} has no motion"


def test_no_css_animation_sits_on_an_ancestor_of_the_fixed_background():
    """A CSS animation creates a containing block, which breaks position:fixed descendants.

    This is why the background field was invisible: `animation: rise` was applied to
    `.block-container > div > div > div > div`, an ancestor of the injected `.bg` element, so
    `position: fixed` silently became `position: absolute` inside a small animated div and the
    field was clipped out of sight. The entrance animation is now scoped to leaf elements.
    """
    css = theme.CSS
    assert ".block-container > div > div > div > div {" not in css, (
        "the entrance animation is back on a generic ancestor chain; anything position:fixed "
        "inside it will be clipped")
    for scoped in (".guide { animation: rise", '[data-testid="stMetric"] { animation: rise',
                   ".stPlotlyChart { animation: rise"):
        assert scoped in css, f"expected the entrance animation scoped to {scoped!r}"


def test_the_background_field_is_actually_visible():
    """It was present, correct and at an opacity that rendered it invisible.

    Two layers at 0.16 and 0.12 under a veil whose centre was 0.94 opaque left nothing to see.
    Pinning a floor means a future tidy-up cannot quietly return it to zero.
    """
    import re

    for selector in (".bg-far", ".bg-near"):
        match = re.search(re.escape(selector) + r"\s*\{[^}]*opacity:\s*([0-9.]+)", theme.CSS)
        assert match, f"{selector} has no opacity declared"
        assert float(match.group(1)) >= 0.3, (
            f"{selector} is at {match.group(1)}, which is not visible on this surface")
    # The FIRST gradient stop is the centre of the page and the only one that can erase the
    # field. A greedy match here captured the last stop instead, which is 0.04, so this check
    # passed regardless of how opaque the middle was. Non-greedy, and asserted on the centre.
    veil = re.search(r"\.bg-veil[^}]*?rgba\(16,20,28,([0-9.]+)\)", theme.CSS)
    assert veil, "the .bg-veil rule is missing entirely; a CSS rewrite deleted it once"
    assert float(veil.group(1)) <= 0.88, (
        f"the veil centre is {veil.group(1)}, opaque enough to erase the field behind the text")
    assert float(veil.group(1)) >= 0.60, (
        f"the veil centre is {veil.group(1)}; the field will run unshaded under the body text")


def test_every_view_has_a_narrator_passage():
    """The guide says what to check. This says why it matters, in the author's voice."""
    import importlib

    app = importlib.import_module("app")
    assert set(app.VIEWS) <= set(dash_copy.NARRATOR), (
        f"views with no narration: {set(app.VIEWS) - set(dash_copy.NARRATOR)}")
    for name, text in dash_copy.NARRATOR.items():
        assert len(text.split()) >= 60, f"{name} narration is too thin to be worth reading"


def test_the_animated_charts_carry_frames_and_a_play_control():
    """Three charts animate. Plotly in Streamlit cannot autoplay, so each needs a control.

    A frame set with no updatemenus is an animation nobody can start, which looks identical to a
    static chart and costs the same to build.
    """
    animated = {
        "animated_model_walk": ("three_arm_models.csv", 3),
        "animated_estimator_walk": ("cross_fitting.csv", 2),
        "animated_protocol_distribution": ("three_arm_per_dataset.csv", 3),
    }
    for name, (table, expected) in animated.items():
        fig = getattr(figures, name)(data.table(table))
        assert len(fig.frames) == expected, f"{name} has {len(fig.frames)} frames"
        assert fig.layout.updatemenus, f"{name} animates but has no play control"
        labels = [b.label for b in fig.layout.updatemenus[0].buttons]
        assert "PLAY" in labels, f"{name} has no PLAY button"
        assert fig.layout.sliders, f"{name} has no step slider"
        assert name in dash_copy.PLAIN, f"{name} has no plain-language caption"


def test_the_guiding_annotations_name_the_finding_on_the_chart():
    """Charts were explained only in surrounding prose, so a reader scanning figures saw shapes.

    These annotations are the argument placed on the artwork. Pinned because they are easy to
    lose in a refactor and their absence is invisible: the chart still renders, it just stops
    saying anything.
    """
    contrast = data.table("three_arm_contrast.csv")
    fig = figures.protocol_contribution(contrast)
    texts = " ".join(a.text for a in fig.layout.annotations)
    assert "lowest contribution" in texts and "highest contribution" in texts, (
        "the contribution chart no longer names the inversion it exists to show")

    null = figures.known_null(data.table("cross_fitting.csv"))
    null_texts = " ".join(a.text for a in null.layout.annotations)
    # Matched case-insensitively on the idea, not on one spelling. The label has been through
    # three forms while the chart's meaning never changed, and a test that pins the exact
    # characters fails on wording rather than on substance.
    assert "zero" in null_texts.lower(), (
        "the known-null chart lost its zero reference; the whole chart is an argument about "
        "what an estimator reports when the answer is zero, so the zero must be labelled")
    assert null.layout.shapes, "the known-null chart lost its zero line"
    assert "cannot exist" in null_texts, "the known-null chart lost its callout"

    scatter = figures.apparent_versus_contribution(data.table("three_arm_per_dataset.csv"))
    quad = " ".join(a.text for a in scatter.layout.annotations)
    assert "contributes most" in quad and "contributes least" in quad, (
        "the scatter lost the quadrant labels that state the inverse relation")
    assert len(scatter.layout.shapes) >= 2, "the scatter lost its quadrant shading"


def test_charts_shown_together_start_their_plot_areas_at_the_same_x():
    """A column of charts only reads as a column if their plot areas line up.

    They did not. The shared template handed every figure a 70px left margin while six charts
    needed 115 to 384 for their own labels, so Plotly grew each one independently and the
    result was a ragged left edge down every view. `_left_margin` now computes the need and
    snaps it to a tier, so any two charts either align exactly or differ obviously.

    Grouped by the view that draws them; charts in different tabs of one view are not compared
    because they are never on screen together.
    """
    groups = {
        "Why this exists": [
            ("survey_timeline", "negative_set_survey_per_method.csv"),
            ("variant_ladder", "variant_ladder.csv"),
        ],
        "How negatives are built": [
            ("animated_matching", "match_quality.csv"),
            ("match_gap_bars", "match_quality.csv"),
            ("matching_cost", "cost_of_matching.csv"),
            ("redraw_stability", "negative_draws.csv"),
        ],
        "Model comparison": [
            ("animated_model_walk", "three_arm_models.csv"),
            ("model_spans", "three_arm_models.csv"),
            ("model_by_protocol", "three_arm_models.csv"),
        ],
        "Protocol sensitivity": [
            ("protocol_levels", "three_arm_contrast.csv"),
            ("protocol_contribution", "three_arm_contrast.csv"),
            ("animated_protocol_distribution", "three_arm_per_dataset.csv"),
            ("apparent_versus_contribution", "three_arm_per_dataset.csv"),
        ],
        "Cross-fitting": [
            ("known_null", "cross_fitting.csv"),
            ("animated_estimator_walk", "cross_fitting.csv"),
            ("estimator_comparison", "cross_fitting.csv"),
        ],
    }
    for view, items in groups.items():
        margins = {}
        for name, table in items:
            fig = getattr(figures, name)(data.table(table))
            margins[name] = fig.layout.margin.l
        assert len(set(margins.values())) == 1, (
            f"{view}: plot areas start at different x positions {margins}")


def test_figures_with_play_controls_reserve_room_above_the_plot():
    """Play buttons sit at paper y above 1. A 46px top margin puts them through the title."""
    for name, table in (("animated_model_walk", "three_arm_models.csv"),
                        ("animated_estimator_walk", "cross_fitting.csv"),
                        ("animated_matching", "match_quality.csv"),
                        ("animated_protocol_distribution", "three_arm_per_dataset.csv")):
        fig = getattr(figures, name)(data.table(table))
        assert fig.layout.margin.t >= 90, f"{name} will clip its play control"
        assert fig.layout.margin.b >= 90, f"{name} will clip its step slider"


def test_the_heatmap_reserves_room_for_its_colourbar():
    """The colourbar is drawn outside the plot area and needs more than the shared right margin."""
    fig = figures.dataset_heatmap(data.table("three_arm_per_dataset.csv"))
    assert fig.layout.margin.r >= 90, "the heatmap colourbar will overlap the plot"
    assert fig.layout.yaxis.showticklabels is False
    assert fig.layout.margin.l == 72, (
        "hidden tick labels should not reserve left margin; 94 dataset names were doing so")


def test_the_declared_filter_state_matches_what_each_view_actually_does():
    """The chip on every view says whether the sidebar filter changes anything. It must be true.

    Reproducibility declared "filters apply in part" while calling `narrow()` zero times, so the
    sidebar did nothing there and the page said otherwise. It was mislabelled because that view
    has its own provenance-class control, which is a different thing from the sidebar filter the
    chip describes.

    Derived from the syntax tree rather than maintained by hand, because this is a claim about
    behaviour and the whole point of the chip is that a reader trusts it.
    """
    import ast

    tree = ast.parse((PACKAGE.parent / "app.py").read_text())
    mismatches = []
    for node in ast.walk(tree):
        if not (isinstance(node, ast.FunctionDef) and node.name.startswith("view_")):
            continue
        narrows = sum(
            1 for child in ast.walk(node)
            if isinstance(child, ast.Call) and getattr(child.func, "id", None) == "narrow"
        )
        title = None
        for child in ast.walk(node):
            if isinstance(child, ast.Call) and getattr(child.func, "id", None) == "page":
                if child.args:
                    title = child.args[0].value
                for keyword in child.keywords:
                    if keyword.arg == "guide":
                        title = keyword.value.value
        if title is None:
            continue
        declared = dash_copy.GUIDE[title]["filters"]
        if narrows == 0 and declared != "static":
            mismatches.append(f"{title!r} says {declared!r} but never calls narrow()")
        if narrows > 0 and declared == "static":
            mismatches.append(f"{title!r} says 'static' but calls narrow() {narrows}x")
    assert not mismatches, "the filter chip lies on:\n  " + "\n  ".join(mismatches)


def test_every_animation_can_be_turned_off():
    """Motion must be suppressible. Two animations had escaped the reduced-motion block.

    `.key-v` was added with an entrance animation and never listed, and before that the block
    named a selector that had stopped animating at all, so it was suppressing nothing. Compared
    on whole selectors here rather than on the last token, because `.rule i` truncated to `i`
    and read as an unsuppressed rule that was in fact covered.
    """
    import re

    css = theme.CSS
    blocks = re.findall(r"@media \(prefers-reduced-motion: reduce\) \{(.*?)\n  \}", css, re.S)
    assert blocks, "there is no reduced-motion block at all"

    suppressed = set()
    for block in blocks:
        for rule in re.finditer(r"([^{}]+)\{([^}]*)\}", block):
            if "animation" in rule.group(2):
                suppressed.update(s.strip() for s in rule.group(1).split(","))

    body = re.sub(r"@media \(prefers-reduced-motion: reduce\) \{.*?\n  \}", "", css, flags=re.S)
    animated = set()
    for rule in re.finditer(r"(?m)^\s*([^{@}\n][^{\n]*?)\s*\{([^}]*)\}", body):
        if re.search(r"animation:\s*(?!none)\w", rule.group(2)):
            animated.update(s.strip() for s in rule.group(1).split(","))

    missing = {a for a in animated if a not in suppressed}
    assert not missing, (
        "these animate but are not disabled under prefers-reduced-motion: " + str(sorted(missing)))


def test_the_background_is_deterministic_and_layered():
    """Four layers at four rates, and the same field on every render.

    The particle layer is placed by seeded arithmetic rather than `random`, because a background
    that reshuffles on each Streamlit rerun is a background the reader notices, which is the one
    thing it must not do.
    """
    first, second = theme.background(), theme.background()
    assert first == second, "the background field changes between renders"
    assert first.count("data:image/svg+xml;base64,") == 4, "expected four layers"
    for layer in ("bg-ring", "bg-dust", "bg-far", "bg-near"):
        assert layer in first, f"{layer} is missing"
    assert 'aria-hidden="true"' in first

    ring = theme._ring(400)
    assert ring.count("<circle") == 3 and ring.count("<line") == 36
    dust = theme._particle_field(600, 400, 20)
    assert dust.count("<text") == 20


def test_no_table_is_registered_without_being_used():
    """A registered table nothing reads is a schema check guarding nothing.

    `class_ratio.csv` sat registered and unread for two releases: the loader validated its
    columns on every run and no view ever drew it, so a regression in it would have failed the
    build while changing nothing a reader sees.
    """
    import ast

    app_src = (PACKAGE.parent / "app.py").read_text()
    fig_src = (PACKAGE / "figures.py").read_text()
    referenced = {
        node.value
        for node in ast.walk(ast.parse(app_src))
        if isinstance(node, ast.Constant) and isinstance(node.value, str)
        and node.value.endswith(".csv")
    }
    unused = sorted(set(data.REQUIRED) - referenced)
    assert not unused, f"registered but never read: {unused}"
    assert "class_ratio" in fig_src or "class_ratio.csv" in app_src


def test_chart_backgrounds_are_set_on_the_layout_not_only_the_template():
    """A consumer that swaps the template takes the backgrounds with it.

    Streamlit's `st.plotly_chart` defaults to `theme="streamlit"`, which REPLACES the figure
    template rather than merging with it. Every chart rendered on white paper inside a dark page
    for several releases because the transparency lived only in the template. Fixed two ways:
    `theme=None` at every call site, and the colours set on the layout as well.
    """
    for name, table in (("protocol_contribution", "three_arm_contrast.csv"),
                        ("known_null", "cross_fitting.csv"),
                        ("dataset_heatmap", "three_arm_per_dataset.csv"),
                        ("protocol_trajectories", "three_arm_per_dataset.csv")):
        fig = getattr(figures, name)(data.table(table))
        assert fig.layout.paper_bgcolor == "rgba(0,0,0,0)", f"{name} has opaque paper"
        assert fig.layout.plot_bgcolor == "rgba(0,0,0,0)", f"{name} has an opaque plot area"

    app_src = (PACKAGE.parent / "app.py").read_text()
    calls = app_src.count("st.plotly_chart(")
    themed = app_src.count("theme=None")
    assert themed >= calls, (
        f"{calls} plotly_chart calls but only {themed} pass theme=None; Streamlit's own theme "
        "will override the figure template on the rest")


def test_the_filter_controls_are_only_drawn_where_they_do_something():
    """Offering a control that does nothing reads as a broken dashboard, not a static view."""
    app_src = (PACKAGE.parent / "app.py").read_text()
    assert "def sidebar_filters(frame, size_column: str, active: bool)" in app_src, (
        "sidebar_filters no longer takes the active flag")
    assert "if not active:" in app_src, "the controls are drawn unconditionally again"
    assert 'uses_filter = copy.GUIDE[GUIDE_KEY.get(choice, choice)]["filters"] != "static"' \
        in app_src, "the sidebar no longer decides from the view's declared filter state"


def test_the_readme_documents_the_palette_that_is_actually_implemented():
    """The README described a palette two revisions old: cyan/amber/rose on #0B0F1A.

    A colour section that names measured separation figures is a claim, and it went stale the
    moment the theme changed. Checked against the module rather than trusted.
    """
    readme = (PACKAGE.parent / "README.md").read_text()
    assert theme.SURFACE in readme, f"README does not name the implemented surface {theme.SURFACE}"
    for colour in theme.CATEGORICAL:
        assert colour in readme, f"README does not name {colour}"
    for stale in ("#0B0F1A", "#22D3EE", "#FB7185"):
        assert stale not in readme, f"README still documents the superseded colour {stale}"


def test_the_readme_does_not_claim_the_page_makes_no_network_request():
    """It loads two font families from Google Fonts, so a blanket claim was false."""
    readme = (PACKAGE.parent / "README.md").read_text()
    assert "reach a network" not in readme, "the blanket no-network claim is back"
    assert "fonts.googleapis.com" in readme, "the font request is no longer disclosed"
    assert "no network request for data" in readme


def test_the_streamlit_config_matches_the_theme_module():
    """`.streamlit/config.toml` colours Streamlit's own widgets. They must be the theme's.

    This drifted: the config kept a superseded neon palette after the module moved to the
    measured one, so sliders and multiselects rendered in the old accent while every chart and
    every surface used the new. Visible only in a browser, which is where it stayed for several
    releases.
    """
    import re

    # REPO ROOT. Streamlit reads $CWD/.streamlit/config.toml and the app is launched from the
    # repository root, so a config under dashboard/ is never read. This test passed against
    # that unread copy, which is exactly how it went unnoticed.
    config_path = data.REPO_ROOT / ".streamlit" / "config.toml"
    assert config_path.exists(), (
        f"{config_path} is missing; a config inside dashboard/ is not read by Streamlit")
    assert not (PACKAGE.parent / ".streamlit" / "config.toml").exists(), (
        "a second config exists under dashboard/, which Streamlit ignores and a reader will "
        "assume is live")
    config = config_path.read_text()

    def value(key):
        match = re.search(rf'^{key}\s*=\s*"([^"]+)"', config, re.M)
        assert match, f"{key} is missing from config.toml"
        return match.group(1)

    assert value("primaryColor") == theme.STEEL
    assert value("backgroundColor") == theme.SURFACE
    assert value("secondaryBackgroundColor") == theme.PANEL
    assert value("textColor") == theme.INK
    assert value("base") == "dark", "the palette was measured against a dark surface"
    assert value("address") == "127.0.0.1", (
        "the server should bind to this machine only; Streamlit otherwise listens on every "
        "interface")


def test_a_view_does_not_repeat_its_own_guide_lead():
    """The explorer printed its lead sentence twice, once in the guide and once below it."""
    app_src = (PACKAGE.parent / "app.py").read_text()
    for name, guide in dash_copy.GUIDE.items():
        lead = guide["lead"]
        # compare on the first clause, which is what a copy-paste duplicate would share
        head = lead.split(",")[0].strip()
        if len(head) < 25:
            continue
        assert app_src.count(head) <= 1, (
            f"{name}: the guide lead {head!r} also appears in the view body")
