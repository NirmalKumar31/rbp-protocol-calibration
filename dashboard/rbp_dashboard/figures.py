"""Plotly figures. Each takes a frame and returns a figure; none reads a file or writes one.

Conventions, applied everywhere rather than per chart: one y-axis, never two. Categorical
colour follows the entity, so removing a protein from a filter does not repaint the survivors.
Error bars are the release's own protein-clustered bootstrap interval, never a recomputed one.
Hover is on by default, because a chart in a browser that cannot be interrogated is a picture.
"""

from __future__ import annotations

import pandas as pd
import plotly.graph_objects as go

from . import theme

TEMPLATE = theme.plotly_template()
ARMS = ("gc", "dn", "neg2")
MODELS = ("kmer", "cnn", "splicebert")


def _base(fig: go.Figure, height: int = 380, ytitle: str = "", xtitle: str = "") -> go.Figure:
    fig.update_layout(template=TEMPLATE, height=height, bargap=0.32, bargroupgap=0.08)
    fig.update_yaxes(title_text=ytitle)
    fig.update_xaxes(title_text=xtitle)
    return fig


def protocol_levels(contrast: pd.DataFrame) -> go.Figure:
    """Composition alone against apparent AUROC, per protocol, with intervals.

    Grouped bars rather than a slope chart: the comparison the reader needs is within a
    protocol, and the two quantities share a scale.
    """
    rows = {
        "Composition alone": "composition alone, {a} arm",
        "Apparent AUROC": "apparent AUROC (composition + score), {a} arm",
    }
    fig = go.Figure()
    for colour, (name, pattern) in zip((theme.SLATE, theme.BLUE), rows.items()):
        vals, los, his = [], [], []
        for arm in ARMS:
            row = contrast[contrast["check"] == pattern.format(a=arm)].iloc[0]
            vals.append(row["value"])
            los.append(row["value"] - row["ci_low"])
            his.append(row["ci_high"] - row["value"])
        fig.add_bar(
            name=name, x=[theme.PROTOCOL_LABEL[a] for a in ARMS], y=vals,
            marker_color=colour, marker_line_width=2, marker_line_color=theme.SURFACE,
            error_y={"type": "data", "array": his, "arrayminus": los,
                     "color": theme.INK_MUTED, "thickness": 1.2, "width": 4},
            hovertemplate="%{x}<br>" + name + " %{y:.3f}<extra></extra>",
        )
    fig.update_yaxes(range=[0.5, 0.92])
    return _base(fig, 390, "AUROC")


def protocol_contribution(contrast: pd.DataFrame) -> go.Figure:
    """The contribution beyond composition, per protocol. Its own chart, its own scale.

    Deliberately not plotted beside the AUROC levels above. They differ by an order of
    magnitude, and a second y-axis is the one chart form this project will not ship.
    """
    vals, los, his, colours = [], [], [], []
    for arm in ARMS:
        row = contrast[contrast["check"] == f"nested contribution, {arm} arm"].iloc[0]
        vals.append(row["value"])
        los.append(row["value"] - row["ci_low"])
        his.append(row["ci_high"] - row["value"])
        colours.append(theme.PROTOCOL_COLOUR[arm])
    fig = go.Figure()
    fig.add_bar(
        x=[theme.PROTOCOL_LABEL[a] for a in ARMS], y=vals,
        marker_color=colours, marker_line_width=2, marker_line_color=theme.SURFACE,
        error_y={"type": "data", "array": his, "arrayminus": los,
                 "color": theme.INK_MUTED, "thickness": 1.2, "width": 5},
        text=[f"{v:.4f}" for v in vals], textposition="outside",
        textfont={"size": 12, "color": theme.INK_MUTED},
        hovertemplate="%{x}<br>contribution %{y:.4f}<extra></extra>",
        showlegend=False,
    )
    return _base(fig, 370, "Nested contribution (AUROC points)")


def contribution_distribution(per_dataset: pd.DataFrame) -> go.Figure:
    """Per-dataset contribution by protocol: every dataset as a point, box behind it.

    The panel means in the bar chart above hide that these distributions overlap. Showing the
    points is the honest version of the same comparison.
    """
    fig = go.Figure()
    for arm in ARMS:
        series = per_dataset[f"gain_{arm}"].dropna()
        fig.add_box(
            name=theme.PROTOCOL_LABEL[arm], y=series,
            marker_color=theme.PROTOCOL_COLOUR[arm], line_width=1.4,
            fillcolor="rgba(0,0,0,0)", boxpoints="all", jitter=0.45, pointpos=0,
            marker={"size": 4, "opacity": 0.55},
            customdata=per_dataset.loc[series.index, "dataset"],
            hovertemplate="%{customdata}<br>contribution %{y:.4f}<extra></extra>",
        )
    fig.add_hline(y=0, line_color=theme.GRID, line_width=1)
    return _base(fig, 420, "Nested contribution (AUROC points)")


def paired_slopes(per_dataset: pd.DataFrame, left: str, right: str, limit: int = 94) -> go.Figure:
    """Dataset-level paired change between two protocols, one thin line per dataset.

    Paired because the datasets are the same in both arms; two independent box plots would
    discard that pairing and overstate the uncertainty of the difference.
    """
    frame = per_dataset.head(limit)
    fig = go.Figure()
    for _, row in frame.iterrows():
        a, b = row[f"gain_{left}"], row[f"gain_{right}"]
        if pd.isna(a) or pd.isna(b):
            continue
        fig.add_scatter(
            x=[theme.PROTOCOL_LABEL[left], theme.PROTOCOL_LABEL[right]], y=[a, b],
            mode="lines", line={"width": 1, "color": theme.BLUE if b > a else theme.AMBER},
            opacity=0.42, showlegend=False, hoverinfo="skip",
        )
    for arm in (left, right):
        vals = frame[f"gain_{arm}"]
        fig.add_scatter(
            x=[theme.PROTOCOL_LABEL[arm]] * len(vals), y=vals, mode="markers",
            marker={"size": 6, "color": theme.PROTOCOL_COLOUR[arm],
                    "line": {"width": 2, "color": theme.SURFACE}},
            customdata=frame["dataset"], showlegend=False,
            hovertemplate="%{customdata}<br>" + theme.PROTOCOL_LABEL[arm]
            + " %{y:.4f}<extra></extra>",
        )
    fig.add_hline(y=0, line_color=theme.GRID, line_width=1)
    return _base(fig, 430, "Nested contribution (AUROC points)")


def apparent_versus_contribution(per_dataset: pd.DataFrame) -> go.Figure:
    """Apparent AUROC against contribution, one series per protocol.

    This is the chart the study's counter-intuitive result lives in: the protocol whose apparent
    AUROC is highest is the one where the model contributes least.
    """
    fig = go.Figure()
    for arm in ARMS:
        fig.add_scatter(
            x=per_dataset[f"full_{arm}"], y=per_dataset[f"gain_{arm}"], mode="markers",
            name=theme.PROTOCOL_LABEL[arm],
            marker={"size": 7, "color": theme.PROTOCOL_COLOUR[arm], "opacity": 0.78,
                    "line": {"width": 2, "color": theme.SURFACE}},
            customdata=per_dataset["dataset"],
            hovertemplate="%{customdata}<br>apparent %{x:.3f}<br>contribution %{y:.4f}"
            "<extra>" + theme.PROTOCOL_LABEL[arm] + "</extra>",
        )
    return _base(fig, 440, "Nested contribution (AUROC points)", "Apparent AUROC")


def model_spans(models: pd.DataFrame) -> go.Figure:
    """Three-protocol span per model class, with intervals. Horizontal, because the labels are long.

    A span of 1.0 is the null, so it is drawn. Without it the bars imply a zero baseline and the
    chart reads as 'CNN is twice the k-mer', which is not what a ratio-to-smallest means.
    """
    fig = go.Figure()
    for model in MODELS:
        row = models[models["check"] == f"{model} three-protocol span"].iloc[0]
        fig.add_bar(
            y=[theme.MODEL_LABEL[model]], x=[row["value"]], orientation="h",
            name=theme.MODEL_LABEL[model], marker_color=theme.MODEL_COLOUR[model],
            marker_line_width=2, marker_line_color=theme.SURFACE,
            error_x={"type": "data", "array": [row["ci_high"] - row["value"]],
                     "arrayminus": [row["value"] - row["ci_low"]],
                     "color": theme.INK_MUTED, "thickness": 1.2, "width": 5},
            text=[f"{row['value']:.2f}x"], textposition="outside",
            textfont={"size": 12, "color": theme.INK_MUTED},
            hovertemplate=theme.MODEL_LABEL[model] + "<br>span %{x:.2f}x<extra></extra>",
            showlegend=False,
        )
    fig.add_vline(x=1.0, line_color=theme.INK_MUTED, line_width=1, line_dash="dot",
                  annotation_text="no protocol effect", annotation_position="top",
                  annotation_font={"size": 11, "color": theme.INK_MUTED})
    fig.update_xaxes(range=[0, 10.4])
    return _base(fig, 330, "", "Span across the three protocols (ratio, largest / smallest)")


def model_by_protocol(models: pd.DataFrame) -> go.Figure:
    """Contribution for each model class within each protocol. Grouped by protocol.

    Grouped this way, not the other, because the claim is about protocols moving the measurement:
    all three model classes put their lowest contribution in the bias-aware arm.
    """
    fig = go.Figure()
    for model in MODELS:
        vals, los, his = [], [], []
        for arm in ARMS:
            row = models[models["check"] == f"{model} nested contribution, {arm} arm"].iloc[0]
            vals.append(row["value"])
            los.append(row["value"] - row["ci_low"])
            his.append(row["ci_high"] - row["value"])
        fig.add_bar(
            name=theme.MODEL_LABEL[model], x=[theme.PROTOCOL_LABEL[a] for a in ARMS], y=vals,
            marker_color=theme.MODEL_COLOUR[model], marker_line_width=2,
            marker_line_color=theme.SURFACE,
            error_y={"type": "data", "array": his, "arrayminus": los,
                     "color": theme.INK_MUTED, "thickness": 1.2, "width": 4},
            hovertemplate="%{x}<br>contribution %{y:.4f}"
            "<extra>" + theme.MODEL_LABEL[model] + "</extra>",
        )
    return _base(fig, 400, "Nested contribution (AUROC points)")


def known_null(cross: pd.DataFrame) -> go.Figure:
    """The 2-mer known-answer test. True contribution is zero by construction.

    The composition baseline already contains every 2-mer frequency, so a 2-mer model can add
    nothing. Whatever an estimator reports here is the estimator measuring itself.
    """
    fig = go.Figure()
    for label, key, colour in (
        ("Two-stage (conventional)", "as published", theme.AMBER),
        ("Cross-fitted (primary)", "fully cross-fitted", theme.BLUE),
    ):
        vals, los, his = [], [], []
        for arm in ARMS:
            row = cross[cross["check"] == f"2-mer contribution {key}, {arm} arm"].iloc[0]
            vals.append(row["value"])
            los.append(max(0.0, row["value"] - row["ci_low"]))
            his.append(max(0.0, row["ci_high"] - row["value"]))
        fig.add_bar(
            name=label, x=[theme.PROTOCOL_LABEL[a] for a in ARMS], y=vals,
            marker_color=colour, marker_line_width=2, marker_line_color=theme.SURFACE,
            error_y={"type": "data", "array": his, "arrayminus": los,
                     "color": theme.INK_MUTED, "thickness": 1.2, "width": 4},
            hovertemplate="%{x}<br>" + label + " %{y:.5f}<extra></extra>",
        )
    fig.add_hline(y=0, line_color=theme.INK_MUTED, line_width=1.4,
                  annotation_text="true value, zero by construction",
                  annotation_position="bottom right",
                  annotation_font={"size": 11, "color": theme.INK_MUTED})
    return _base(fig, 400, "Reported contribution (AUROC points)")


def estimator_comparison(cross: pd.DataFrame) -> go.Figure:
    """The 4-mer under both estimators, per protocol. The real signal, for contrast with the null.

    Unlike the 2-mer, cross-fitting does not collapse these. It moves them slightly upward,
    which is why the paper reports the cross-fitted value as primary rather than as a correction.
    """
    fig = go.Figure()
    for label, key, colour in (
        ("Two-stage (conventional)", "as published", theme.AMBER),
        ("Cross-fitted (primary)", "fully cross-fitted", theme.BLUE),
    ):
        vals, los, his = [], [], []
        for arm in ARMS:
            row = cross[cross["check"] == f"4-mer contribution {key}, {arm} arm"].iloc[0]
            vals.append(row["value"])
            los.append(row["value"] - row["ci_low"])
            his.append(row["ci_high"] - row["value"])
        fig.add_bar(
            name=label, x=[theme.PROTOCOL_LABEL[a] for a in ARMS], y=vals,
            marker_color=colour, marker_line_width=2, marker_line_color=theme.SURFACE,
            error_y={"type": "data", "array": his, "arrayminus": los,
                     "color": theme.INK_MUTED, "thickness": 1.2, "width": 4},
            hovertemplate="%{x}<br>" + label + " %{y:.4f}<extra></extra>",
        )
    return _base(fig, 400, "4-mer nested contribution (AUROC points)")


def external_arms(per_dataset: pd.DataFrame) -> go.Figure:
    """The held-out benchmark's two negative sets, per dataset.

    Their labels are negative-1 and negative-2, not this study's three protocols, and are not
    renamed here. Mapping them onto GC-matched or bias-aware would assert an equivalence the
    deposit does not support.
    """
    fig = go.Figure()
    for arm, label, colour in (
        ("n1", "negative-1 (bias-agnostic)", theme.BLUE),
        ("n2", "negative-2 (bias-aware)", theme.AMBER),
    ):
        series = per_dataset[f"gain_{arm}"].dropna()
        fig.add_box(
            name=label, y=series, marker_color=colour, line_width=1.4,
            fillcolor="rgba(0,0,0,0)", boxpoints="all", jitter=0.45, pointpos=0,
            marker={"size": 4, "opacity": 0.5},
            customdata=per_dataset.loc[series.index, "dataset"],
            hovertemplate="%{customdata}<br>contribution %{y:.4f}<extra></extra>",
        )
    fig.add_hline(y=0, line_color=theme.GRID, line_width=1)
    return _base(fig, 400, "Nested contribution (AUROC points)")


def external_levels(external: pd.DataFrame, which: str) -> go.Figure:
    """One measure, one chart, one axis. `which` is "baseline" or "contribution".

    Drawn as two charts rather than one because the two measures differ by an order of
    magnitude. Putting them together needs either a second y-axis or a scaling factor on one
    series, and a scaling factor is a second axis wearing a disguise.

    This is the finding that did NOT replicate. Here the arm with the higher composition
    baseline also has the higher contribution, so the two move together. In the 94-dataset panel
    they move in opposite directions.
    """
    labels = {
        "baseline": (
            "composition baseline, negative-{i}",
            "Composition baseline (AUROC)", theme.SLATE, "{v:.3f}", 3,
        ),
        "contribution": (
            "nested contribution, negative-{i} ({tag})",
            "Nested contribution (AUROC points)", theme.BLUE, "{v:.4f}", 4,
        ),
    }
    pattern, ytitle, colour, fmt, digits = labels[which]
    tags = {"1": "bias-agnostic", "2": "bias-aware"}
    names, vals, los, his = [], [], [], []
    for i in ("1", "2"):
        row = external[external["check"] == pattern.format(i=i, tag=tags[i])].iloc[0]
        names.append(f"negative-{i}\n({tags[i]})")
        vals.append(row["value"])
        los.append(row["value"] - row["ci_low"])
        his.append(row["ci_high"] - row["value"])
    fig = go.Figure()
    fig.add_bar(
        x=names, y=vals, marker_color=colour, marker_line_width=2,
        marker_line_color=theme.SURFACE,
        error_y={"type": "data", "array": his, "arrayminus": los,
                 "color": theme.INK_MUTED, "thickness": 1.2, "width": 5},
        text=[fmt.format(v=v) for v in vals], textposition="outside",
        textfont={"size": 12, "color": theme.INK_MUTED},
        hovertemplate="%{x}<br>" + f"%{{y:.{digits}f}}" + "<extra></extra>",
        showlegend=False,
    )
    fig.update_yaxes(range=[0, max(vals) * 1.3])
    return _base(fig, 360, ytitle)


def provenance_mix(provenance: pd.DataFrame) -> go.Figure:
    """How the committed tables were produced. Sorted by count, one bar per class.

    A single colour: these are categories of one measure, not competing series, and colouring
    them differently would imply a ranking the classes do not carry.
    """
    counts = provenance["status"].value_counts().sort_values()
    fig = go.Figure()
    fig.add_bar(
        y=counts.index.tolist(), x=counts.to_list(), orientation="h",
        marker_color=theme.BLUE, marker_line_width=2, marker_line_color=theme.SURFACE,
        text=counts.to_list(), textposition="outside",
        textfont={"size": 12, "color": theme.INK_MUTED},
        hovertemplate="%{y}<br>%{x} tables<extra></extra>", showlegend=False,
    )
    fig.update_xaxes(range=[0, max(counts) * 1.18])
    return _base(fig, 330, "", "Committed tables")
