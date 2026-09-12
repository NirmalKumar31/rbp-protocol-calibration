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


# One column of plot areas down the page only looks like one column if every figure reserves the
# same left margin. This is the floor; a figure whose labels need more gets more, and the audit
# that produced these numbers found six charts asking for 115 to 384 pixels while the template
# handed every one of them 70.
_MARGIN_L_MIN = 72
_MARGIN_L_MAX = 260
_MARGIN_TIERS = (72, 136, 208, 260)
_CHAR_PX = 6.35          # measured against JetBrains Mono at the tick size used here
_LABEL_PAD = 24


def _left_margin(fig: go.Figure) -> int:
    """Widest categorical y label, converted to pixels.

    Plotly does grow the margin on its own, but it does so per figure and only far enough for
    that figure, which is precisely what makes a column of charts look ragged. Computing it here
    and clamping to a shared floor keeps the plot areas aligned with each other instead.
    """
    if fig.layout.yaxis.showticklabels is False:
        return _MARGIN_L_MIN
    labels: list[str] = []
    for trace in fig.data:
        ticks = getattr(trace, "y", None)
        if ticks is not None:
            labels += [str(v) for v in ticks if isinstance(v, str)]
    ticktext = fig.layout.yaxis.ticktext
    if ticktext:
        labels += [str(v) for v in ticktext]
    if not labels:
        return _MARGIN_L_MIN
    widest = max(len(label) for label in labels)
    needed = widest * _CHAR_PX + _LABEL_PAD
    # Snap to tiers rather than using the exact width. Five charts came out at 157, 182, 189,
    # 189 and 201, which is close enough to read as a mistake rather than as a decision: on a
    # page with two of them the plot areas start nineteen pixels apart. Tiers mean any two
    # charts either align exactly or differ obviously.
    for tier in _MARGIN_TIERS:
        if needed <= tier:
            return tier
    return _MARGIN_L_MAX


def _base(fig: go.Figure, height: int = 380, ytitle: str = "", xtitle: str = "") -> go.Figure:
    """Apply the shared template, then reserve the space this particular figure needs.

    Top margin is decided by whether the figure carries play controls: those sit above the plot
    at paper y above 1, so a 40px top clips them into the title.
    """
    has_controls = bool(fig.layout.updatemenus)
    has_slider = bool(fig.layout.sliders)
    fig.update_layout(
        template=TEMPLATE, height=height, bargap=0.32, bargroupgap=0.08,
        # Set on the layout as well as in the template. A consumer that swaps the template
        # takes the backgrounds with it: Streamlit's default `theme="streamlit"` did exactly
        # that and rendered every chart on white paper inside a dark page.
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        margin={
            "l": _left_margin(fig),
            "r": 36,
            "t": 96 if has_controls else 46,
            "b": 92 if has_slider else 58,
        },
    )
    fig.update_yaxes(title_text=ytitle, automargin=True)
    fig.update_xaxes(title_text=xtitle, automargin=True)
    return fig


def _arrow(fig, x, y, text: str, ax: int = 40, ay: int = -46, colour: str | None = None):
    """An annotation that points. Used to name the finding on the chart itself.

    Charts were being explained entirely in surrounding prose, which means a reader scanning the
    figures sees shapes and no argument. One arrow per chart, at most, on the feature the view
    is actually about.
    """
    fig.add_annotation(
        x=x, y=y, text=text, showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=1.1,
        arrowcolor=colour or theme.INK_FAINT, ax=ax, ay=ay,
        font={"size": 10.5, "color": colour or theme.INK_MUTED, "family": theme.FONT},
        align="left", bgcolor="rgba(21,26,36,0.92)", bordercolor=theme.BORDER,
        borderwidth=1, borderpad=5,
    )
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
    for colour, (name, pattern) in zip((theme.SLATE, theme.STEEL), rows.items()):
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
    fig.update_yaxes(range=[0.5, 0.95])
    fig.add_annotation(
        x=theme.PROTOCOL_LABEL["dn"], y=0.72, showarrow=False,
        text="the gap is<br>what the model adds", xanchor="center",
        font={"size": 10, "color": theme.INK_FAINT, "family": theme.FONT},
        bgcolor="rgba(21,26,36,0.9)", bordercolor=theme.BORDER, borderwidth=1, borderpad=4)
    return _base(fig, 410, "AUROC")


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
    top = max(vals)
    _arrow(fig, theme.PROTOCOL_LABEL["neg2"], vals[2],
           "lowest contribution,<br>highest headline score", ax=-4, ay=-72, colour=theme.SLATE)
    _arrow(fig, theme.PROTOCOL_LABEL["dn"], top,
           "highest contribution,<br>lowest headline score", ax=10, ay=-38, colour=theme.BRASS)
    fig.update_yaxes(range=[0, top * 1.55])
    return _base(fig, 400, "Nested contribution (AUROC points)")


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
            mode="lines", line={"width": 1, "color": theme.STEEL if b > a else theme.BRASS},
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
    # The shading is the claim: high apparent score sits with low contribution and the reverse.
    fig.add_shape(type="rect", xref="paper", yref="paper", x0=0, x1=0.42, y0=0.42, y1=1,
                  fillcolor=theme.BRASS, opacity=0.05, line_width=0, layer="below")
    fig.add_shape(type="rect", xref="paper", yref="paper", x0=0.58, x1=1, y0=0, y1=0.5,
                  fillcolor=theme.SLATE, opacity=0.06, line_width=0, layer="below")
    fig.add_annotation(
        xref="paper", yref="paper", x=0.02, y=0.96, showarrow=False,
        text="looks hard, contributes most", xanchor="left",
        font={"size": 10, "color": theme.BRASS, "family": theme.MONO})
    fig.add_annotation(
        xref="paper", yref="paper", x=0.98, y=0.04, showarrow=False,
        text="looks easy, contributes least", xanchor="right",
        font={"size": 10, "color": theme.SLATE, "family": theme.MONO})
    return _base(fig, 470, "Nested contribution (AUROC points)", "Apparent AUROC")


def model_spans(models: pd.DataFrame) -> go.Figure:
    """Three-protocol span per model class, with intervals. Horizontal, because the labels are long.

    A span of 1.0 is the null, so it is drawn. Without it the bars imply a zero baseline and the
    chart reads as 'CNN is twice the k-mer', which is not what a ratio-to-smallest means.
    """
    fig = go.Figure()
    for model in MODELS:
        row = models[models["check"] == f"{model} three-protocol span"].iloc[0]
        # Short tick labels, full name on hover. The long form pushed this chart's plot area
        # 136px right of the two vertical charts beside it on the same view.
        fig.add_bar(
            y=[theme.MODEL_SHORT[model]], x=[row["value"]], orientation="h",
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
    # Names inside the plot rather than on the axis. With three bars the axis labels bought
    # nothing and pushed this chart's plot area 64px right of the two vertical charts sharing
    # its view, which is the sort of misalignment that reads as carelessness.
    for model in MODELS:
        fig.add_annotation(
            x=0.12, y=theme.MODEL_SHORT[model], text=theme.MODEL_LABEL[model],
            showarrow=False, xanchor="left", yanchor="bottom", yshift=13,
            font={"size": 10.5, "color": theme.INK_MUTED, "family": theme.FONT},
        )
    fig.update_yaxes(showticklabels=False)
    fig.update_xaxes(range=[0, 10.4])
    return _base(fig, 360, "", "Span across the three protocols (ratio, largest / smallest)")


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
        ("Two-stage (conventional)", "as published", theme.BRASS),
        ("Cross-fitted (primary)", "fully cross-fitted", theme.STEEL),
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
    worst = max(
        cross[cross["check"] == f"2-mer contribution as published, {a} arm"].iloc[0]["value"]
        for a in ARMS
    )
    # A line, not a band. The band was two percent of the axis tall and read as a third
    # series sitting on the floor, and the label inside the plot collided with the first pair
    # of bars. Both now live outside the data area.
    # Labelled on the line itself, at the right end where no bar reaches. Anchoring it outside
    # the plot needed a wider left margin, which pushed this chart's plot area out of line with
    # the two beside it on the same view.
    fig.add_hline(
        y=0, line_color=theme.INK, line_width=1.4, opacity=0.65,
        annotation_text="true value, zero by construction",
        annotation_position="top right",
        annotation_font={"size": 10, "color": theme.INK_MUTED, "family": theme.MONO},
    )
    _arrow(fig, theme.PROTOCOL_LABEL["dn"], worst,
           "the conventional estimator<br>reporting signal that cannot exist",
           ax=26, ay=-44, colour=theme.BRASS)
    fig.update_yaxes(range=[-worst * 0.34, worst * 1.5])
    return _base(fig, 430, "Reported contribution (AUROC points)")


def estimator_comparison(cross: pd.DataFrame) -> go.Figure:
    """The 4-mer under both estimators, per protocol. The real signal, for contrast with the null.

    Unlike the 2-mer, cross-fitting does not collapse these. It moves them slightly upward,
    which is why the paper reports the cross-fitted value as primary rather than as a correction.
    """
    fig = go.Figure()
    for label, key, colour in (
        ("Two-stage (conventional)", "as published", theme.BRASS),
        ("Cross-fitted (primary)", "fully cross-fitted", theme.STEEL),
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
        ("n1", "negative-1 (bias-agnostic)", theme.STEEL),
        ("n2", "negative-2 (bias-aware)", theme.BRASS),
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
            "Nested contribution (AUROC points)", theme.STEEL, "{v:.4f}", 4,
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
        marker_color=theme.STEEL, marker_line_width=2, marker_line_color=theme.SURFACE,
        text=counts.to_list(), textposition="outside",
        textfont={"size": 12, "color": theme.INK_MUTED},
        hovertemplate="%{y}<br>%{x} tables<extra></extra>", showlegend=False,
    )
    fig.update_xaxes(range=[0, max(counts) * 1.18])
    return _base(fig, 330, "", "Committed tables")


def dataset_heatmap(per_dataset: pd.DataFrame, sort_by: str = "dn") -> go.Figure:
    """Every dataset against the three protocols, as a sequential heatmap.

    One hue, light to dark, because this is magnitude and not identity. A diverging scale would
    imply a meaningful midpoint and a rainbow would imply an order the values do not carry.
    Sorted so the structure is visible rather than alphabetical.
    """
    frame = per_dataset.sort_values(f"gain_{sort_by}", ascending=True)
    z = [[row[f"gain_{a}"] for a in ARMS] for _, row in frame.iterrows()]
    fig = go.Figure(
        go.Heatmap(
            z=z, x=[theme.PROTOCOL_SHORT[a] for a in ARMS], y=frame["dataset"].tolist(),
            colorscale=[[0, "#12171F"], [0.3, "#2B3A48"], [0.62, "#5B7A93"], [1, "#A8C8E0"]],
            hovertemplate="%{y}<br>%{x}  %{z:.4f}<extra></extra>",
            colorbar={
                "title": {"text": "contribution", "font": {"size": 11, "color": theme.INK_MUTED}},
                "tickfont": {"size": 10, "color": theme.INK_MUTED},
                "outlinewidth": 0, "thickness": 12, "len": 0.72,
            },
        )
    )
    fig.update_yaxes(showticklabels=False, title_text=f"{len(frame)} datasets")
    fig = _base(fig, 560, "", "")
    # The colourbar is drawn outside the plot area, so the shared right margin is not enough.
    fig.update_layout(margin_r=96)
    return fig


def ranked_datasets(per_dataset: pd.DataFrame, arm: str, top: int = 24) -> go.Figure:
    """The strongest datasets in one protocol. A bar chart because the job is ranking.

    Truncated rather than scrolled: 94 labelled bars is a wall, and the tail is visible in the
    heatmap and the distribution already.
    """
    frame = per_dataset.nlargest(top, f"gain_{arm}")
    fig = go.Figure()
    fig.add_bar(
        y=frame["dataset"], x=frame[f"gain_{arm}"], orientation="h",
        marker_color=theme.PROTOCOL_COLOUR[arm], marker_line_width=2,
        marker_line_color=theme.SURFACE,
        customdata=frame[["protein", "cell"]].to_numpy(),
        hovertemplate="%{y}<br>%{customdata[0]} in %{customdata[1]}"
                      "<br>contribution %{x:.4f}<extra></extra>",
        showlegend=False,
    )
    fig.update_yaxes(autorange="reversed", tickfont={"size": 10})
    return _base(fig, 34 * min(top, len(frame)) + 90, "",
                 f"Nested contribution, {theme.PROTOCOL_LABEL[arm]}")


def dataset_profile(models_row: pd.Series) -> go.Figure:
    """One dataset: every model class against every protocol. Grouped bars, one axis.

    A radar chart would be the obvious choice and the wrong one: three axes with a shared scale
    invite reading the enclosed area as a quantity, and the area of a radar polygon means
    nothing.
    """
    fig = go.Figure()
    for model in MODELS:
        vals = [models_row.get(f"{model}_gain_{a}") for a in ARMS]
        fig.add_bar(
            name=theme.MODEL_LABEL[model], x=[theme.PROTOCOL_SHORT[a] for a in ARMS], y=vals,
            marker_color=theme.MODEL_COLOUR[model], marker_line_width=2,
            marker_line_color=theme.SURFACE,
            text=[f"{v:.3f}" if v is not None else "" for v in vals],
            textposition="outside", textfont={"size": 10, "color": theme.INK_MUTED},
            hovertemplate="%{x}<br>%{y:.4f}<extra>" + theme.MODEL_LABEL[model] + "</extra>",
        )
    fig.add_hline(y=0, line_color=theme.BORDER, line_width=1)
    return _base(fig, 360, "Nested contribution (AUROC points)")


def animated_model_walk(models: pd.DataFrame) -> go.Figure:
    """Plotly frames stepping through the three model classes. A real animation, not a CSS fade.

    Worth the machinery for one chart: the claim is that the protocol ordering survives a change
    of model class, and watching the bars move while the ordering holds carries that better than
    three static panels do.
    """
    frames, buttons = [], []
    for model in MODELS:
        vals = [models[models["check"] == f"{model} nested contribution, {a} arm"].iloc[0]["value"]
                for a in ARMS]
        frames.append(go.Frame(
            name=theme.MODEL_SHORT[model],
            data=[go.Bar(
                x=[theme.PROTOCOL_LABEL[a] for a in ARMS], y=vals,
                marker_color=[theme.PROTOCOL_COLOUR[a] for a in ARMS],
                marker_line_width=2, marker_line_color=theme.SURFACE,
                text=[f"{v:.4f}" for v in vals], textposition="outside",
                textfont={"size": 11, "color": theme.INK_MUTED},
                hovertemplate="%{x}<br>%{y:.4f}<extra></extra>",
            )],
        ))
        buttons.append({
            "label": theme.MODEL_SHORT[model], "method": "animate",
            "args": [[theme.MODEL_SHORT[model]],
                     {"mode": "immediate", "frame": {"duration": 500, "redraw": True},
                      "transition": {"duration": 450, "easing": "cubic-in-out"}}],
        })

    fig = go.Figure(data=frames[0].data, frames=frames)
    fig.update_layout(
        updatemenus=[{
            "type": "buttons", "direction": "right", "showactive": True,
            "x": 0, "y": 1.22, "xanchor": "left", "yanchor": "top",
            "pad": {"r": 6, "t": 0},
            "bgcolor": theme.RAISED, "bordercolor": theme.BORDER,
            "font": {"family": theme.MONO, "size": 11, "color": theme.INK_MUTED},
            "buttons": [{
                "label": "PLAY", "method": "animate",
                "args": [None, {"frame": {"duration": 900, "redraw": True},
                                "transition": {"duration": 450, "easing": "cubic-in-out"},
                                "fromcurrent": True}],
            }] + buttons,
        }],
        sliders=[{
            "active": 0, "x": 0, "len": 1.0, "y": -0.16,
            "pad": {"t": 28},
            "currentvalue": {"prefix": "model class  ", "font": {
                "family": theme.MONO, "size": 12, "color": theme.STEEL}},
            "tickcolor": theme.BORDER,
            "font": {"family": theme.MONO, "size": 10, "color": theme.INK_MUTED},
            "steps": [{
                "label": theme.MODEL_SHORT[m], "method": "animate",
                "args": [[theme.MODEL_SHORT[m]],
                         {"mode": "immediate", "frame": {"duration": 500, "redraw": True},
                          "transition": {"duration": 450, "easing": "cubic-in-out"}}],
            } for m in MODELS],
        }],
    )
    fig.update_yaxes(range=[0, 0.205])
    return _base(fig, 470, "Nested contribution (AUROC points)")


def filtered_versus_published(
    subset_stat: dict | None, published: tuple[float, float | None, float | None], arm: str
) -> go.Figure:
    """A filtered subset's descriptive mean beside the published panel estimate.

    Drawn together because the comparison is the point, and drawn differently because they are
    different kinds of thing. The published value carries its bootstrap interval; the live one
    carries none, because a filtered mean has no clustered interval and inventing one would be
    the exact defect this project keeps finding.
    """
    val, low, high = published
    fig = go.Figure()
    fig.add_bar(
        name="Published panel estimate", x=["published"], y=[val],
        marker_color=theme.PROTOCOL_COLOUR[arm], marker_line_width=2,
        marker_line_color=theme.SURFACE,
        error_y=({"type": "data", "array": [high - val], "arrayminus": [val - low],
                  "color": theme.INK_MUTED, "thickness": 1.3, "width": 6}
                 if low is not None and high is not None else None),
        text=[f"{val:.4f}"], textposition="outside",
        textfont={"size": 11, "color": theme.INK_MUTED},
        hovertemplate="published %{y:.4f}, 94 datasets<extra></extra>",
    )
    if subset_stat:
        fig.add_bar(
            name="Filtered subset (descriptive)", x=["filtered"], y=[subset_stat["mean"]],
            marker_color=theme.INK_FAINT, marker_line_width=2,
            marker_line_color=theme.SURFACE,
            marker_pattern={"shape": "/", "size": 4, "solidity": 0.22},
            text=[f"{subset_stat['mean']:.4f}"], textposition="outside",
            textfont={"size": 11, "color": theme.INK_MUTED},
            hovertemplate=f"descriptive mean of {subset_stat['n']} datasets"
                          "<br>%{y:.4f}, no interval<extra></extra>",
        )
    return _base(fig, 330, "Nested contribution (AUROC points)")


def _play_controls(labels, y=1.24, prefix="frame  "):
    """A PLAY button and a step slider, styled once.

    Plotly in Streamlit cannot autoplay on load, so every animation here is deliberate: the
    reader presses play. That is a constraint rather than a choice, and it shapes which charts
    get animated at all. Motion is only worth a control where the movement itself is the
    argument, so there are two of these and not ten.
    """
    steps = [{
        "label": name, "method": "animate",
        "args": [[name], {"mode": "immediate", "frame": {"duration": 520, "redraw": True},
                          "transition": {"duration": 460, "easing": "cubic-in-out"}}],
    } for name in labels]
    return (
        [{
            "type": "buttons", "direction": "right", "showactive": False,
            "x": 0, "y": y, "xanchor": "left", "yanchor": "top", "pad": {"r": 6, "t": 0},
            "bgcolor": theme.RAISED, "bordercolor": theme.BORDER,
            "font": {"family": theme.MONO, "size": 10.5, "color": theme.INK_MUTED},
            "buttons": [{
                "label": "PLAY", "method": "animate",
                "args": [None, {"frame": {"duration": 1000, "redraw": True},
                                "transition": {"duration": 460, "easing": "cubic-in-out"},
                                "fromcurrent": True}],
            }] + steps,
        }],
        [{
            "active": 0, "x": 0, "len": 1.0, "y": -0.17, "pad": {"t": 28},
            "currentvalue": {"prefix": prefix, "font": {
                "family": theme.MONO, "size": 11.5, "color": theme.STEEL}},
            "tickcolor": theme.BORDER,
            "font": {"family": theme.MONO, "size": 9.5, "color": theme.INK_FAINT},
            "steps": steps,
        }],
    )


def animated_estimator_walk(cross: pd.DataFrame) -> go.Figure:
    """The span narrowing as the estimator is corrected. Press play.

    The movement is the argument here: the three bars shift, and the ratio between the tallest
    and the shortest goes from 5.42 to 4.84 while their ordering holds. A pair of static panels
    shows the same numbers and not the same thing.
    """
    stages = [
        ("TWO-STAGE", "as published", "4-mer three-arm span, as published", theme.BRASS),
        ("CROSS-FITTED", "fully cross-fitted", "4-mer three-arm span, fully cross-fitted",
         theme.STEEL),
    ]
    frames, ceiling = [], 0.0
    for name, key, span_key, colour in stages:
        vals = [cross[cross["check"] == f"4-mer contribution {key}, {a} arm"].iloc[0]["value"]
                for a in ARMS]
        ceiling = max(ceiling, max(vals))
        span = cross[cross["check"] == span_key].iloc[0]["value"]
        frames.append(go.Frame(
            name=name,
            data=[go.Bar(
                x=[theme.PROTOCOL_LABEL[a] for a in ARMS], y=vals,
                marker_color=colour, marker_line_width=2, marker_line_color=theme.SURFACE,
                text=[f"{v:.4f}" for v in vals], textposition="outside",
                textfont={"size": 11, "color": theme.INK_MUTED},
                hovertemplate="%{x}<br>%{y:.4f}<extra>" + name + "</extra>",
            )],
            layout=go.Layout(annotations=[{
                "xref": "paper", "yref": "paper", "x": 0.99, "y": 0.97, "showarrow": False,
                "text": f"span  {span:.2f}x", "xanchor": "right",
                "font": {"family": theme.MONO, "size": 15, "color": colour},
            }]),
        ))

    fig = go.Figure(data=frames[0].data, frames=frames)
    menus, sliders = _play_controls([f.name for f in frames], prefix="estimator  ")
    fig.update_layout(updatemenus=menus, sliders=sliders,
                      annotations=list(frames[0].layout.annotations))
    fig.update_yaxes(range=[0, ceiling * 1.32])
    return _base(fig, 470, "4-mer nested contribution (AUROC points)")


def animated_protocol_distribution(per_dataset: pd.DataFrame) -> go.Figure:
    """One protocol at a time, as a histogram. Press play to step through the three.

    Overlaying three histograms produces a mess at this sample size, and three small multiples
    lose the shared axis that makes the shift legible. Stepping keeps one axis and one shape on
    screen at a time.
    """
    edges = 26
    lo = min(per_dataset[f"gain_{a}"].min() for a in ARMS)
    hi = max(per_dataset[f"gain_{a}"].max() for a in ARMS)
    frames, tallest = [], 0
    for arm in ARMS:
        series = per_dataset[f"gain_{arm}"].dropna()
        counts, _ = pd.cut(series, bins=edges, retbins=True)
        tallest = max(tallest, counts.value_counts().max())
        frames.append(go.Frame(
            name=theme.PROTOCOL_SHORT[arm],
            data=[go.Histogram(
                x=series, nbinsx=edges, marker_color=theme.PROTOCOL_COLOUR[arm],
                marker_line_width=1.5, marker_line_color=theme.SURFACE, opacity=0.9,
                hovertemplate="%{x:.4f}<br>%{y} datasets<extra></extra>",
            )],
            layout=go.Layout(annotations=[{
                "xref": "paper", "yref": "paper", "x": 0.99, "y": 0.96, "showarrow": False,
                "text": f"{theme.PROTOCOL_LABEL[arm]}<br>median {series.median():.4f}",
                "xanchor": "right", "align": "right",
                "font": {"family": theme.MONO, "size": 11,
                         "color": theme.PROTOCOL_COLOUR[arm]},
            }]),
        ))
    fig = go.Figure(data=frames[0].data, frames=frames)
    menus, sliders = _play_controls([f.name for f in frames], prefix="protocol  ")
    fig.update_layout(updatemenus=menus, sliders=sliders,
                      annotations=list(frames[0].layout.annotations), bargap=0.06)
    fig.update_xaxes(range=[lo * 1.05, hi * 1.05])
    fig.update_yaxes(range=[0, tallest * 1.22])
    return _base(fig, 450, "Datasets", "Nested contribution (AUROC points)")


def match_quality_curve(match: pd.DataFrame) -> go.Figure:
    """How closely each protocol actually matched composition. The construction, measured.

    Three step curves: the share of bound/unbound pairs whose GC fraction differs by less than a
    threshold. This is the mechanism the rest of the study is about. The two composition-matched
    arms rise steeply; the bias-aware arm does not, because it is not matching composition at
    all, and that difference is what moves its baseline and therefore its measured contribution.
    """
    thresholds = [0.0, 0.05, 0.10, 0.15]
    fig = go.Figure()
    for arm in ARMS:
        ys = [0.0] + [
            match[match["check"] == f"fraction of pairs within |dGC| {t:.2f}, {arm} arm"]
            .iloc[0]["value"] for t in thresholds[1:]
        ]
        fig.add_scatter(
            x=thresholds, y=ys, mode="lines+markers", name=theme.PROTOCOL_LABEL[arm],
            line={"width": 2.2, "color": theme.PROTOCOL_COLOUR[arm], "shape": "spline"},
            marker={"size": 8, "color": theme.PROTOCOL_COLOUR[arm],
                    "line": {"width": 2, "color": theme.SURFACE}},
            hovertemplate="within %{x:.2f} GC<br>%{y:.1%} of pairs"
                          "<extra>" + theme.PROTOCOL_LABEL[arm] + "</extra>",
        )
    fig.add_hline(y=1.0, line_color=theme.BORDER, line_width=1, line_dash="dot")
    fig.add_annotation(
        x=0.05, y=0.204, ax=54, ay=34, text="only 20% of bias-aware pairs<br>match on composition",
        showarrow=True, arrowhead=2, arrowwidth=1.1, arrowcolor=theme.SLATE,
        font={"size": 10.5, "color": theme.SLATE}, align="left",
        bgcolor="rgba(21,26,36,0.92)", bordercolor=theme.BORDER, borderwidth=1, borderpad=5)
    fig.update_yaxes(range=[0, 1.08], tickformat=".0%")
    return _base(fig, 420, "Share of pairs matched this closely",
                 "Allowed difference in GC fraction")


def match_gap_bars(match: pd.DataFrame) -> go.Figure:
    """Median composition gap per arm, with the tail behind it. Lower is a tighter match."""
    fig = go.Figure()
    for label, key, opacity in (("median", "gc_gap_median", 1.0),
                                ("90th percentile", "gc_gap_p90", 0.45)):
        vals = [match[match["check"] == f"{key}, {a} arm"].iloc[0]["value"] for a in ARMS]
        fig.add_bar(
            name=label, x=[theme.PROTOCOL_LABEL[a] for a in ARMS], y=vals,
            marker_color=[theme.PROTOCOL_COLOUR[a] for a in ARMS], marker_opacity=opacity,
            marker_line_width=2, marker_line_color=theme.SURFACE,
            text=[f"{v:.3f}" for v in vals], textposition="outside",
            textfont={"size": 10, "color": theme.INK_FAINT},
            hovertemplate="%{x}<br>" + label + " |dGC| %{y:.4f}<extra></extra>",
        )
    return _base(fig, 380, "Difference in GC fraction, bound vs unbound")


def animated_matching(match: pd.DataFrame) -> go.Figure:
    """Step through the three protocols, watching the match tighten or loosen. Press play.

    The same curve, one arm at a time, on a fixed axis. Stepping rather than overlaying because
    the point is the change between them and a reader tracking three splines at once loses it.
    """
    thresholds = [0.0, 0.05, 0.10, 0.15]
    frames = []
    for arm in ARMS:
        ys = [0.0] + [
            match[match["check"] == f"fraction of pairs within |dGC| {t:.2f}, {arm} arm"]
            .iloc[0]["value"] for t in thresholds[1:]
        ]
        med = match[match["check"] == f"gc_gap_median, {arm} arm"].iloc[0]["value"]
        frames.append(go.Frame(
            name=theme.PROTOCOL_SHORT[arm],
            data=[go.Scatter(
                x=thresholds, y=ys, mode="lines+markers", fill="tozeroy",
                line={"width": 2.4, "color": theme.PROTOCOL_COLOUR[arm], "shape": "spline"},
                fillcolor=theme.FILL_TINT[arm],
                marker={"size": 9, "color": theme.PROTOCOL_COLOUR[arm],
                        "line": {"width": 2, "color": theme.SURFACE}},
                hovertemplate="within %{x:.2f} GC<br>%{y:.1%} of pairs<extra></extra>",
            )],
            layout=go.Layout(annotations=[{
                "xref": "paper", "yref": "paper", "x": 0.98, "y": 0.12, "showarrow": False,
                "text": f"{theme.PROTOCOL_LABEL[arm]}<br>median gap {med:.4f}",
                "xanchor": "right", "align": "right",
                "font": {"family": theme.MONO, "size": 11,
                         "color": theme.PROTOCOL_COLOUR[arm]},
            }]),
        ))
    fig = go.Figure(data=frames[0].data, frames=frames)
    menus, sliders = _play_controls([f.name for f in frames], prefix="protocol  ")
    fig.update_layout(updatemenus=menus, sliders=sliders, showlegend=False,
                      annotations=list(frames[0].layout.annotations))
    fig.update_yaxes(range=[0, 1.08], tickformat=".0%")
    fig.update_xaxes(range=[0, 0.155])
    return _base(fig, 450, "Share of pairs matched this closely",
                 "Allowed difference in GC fraction")


def matching_cost(cost: pd.DataFrame) -> go.Figure:
    """What stricter matching costs the model's raw score, per dataset.

    Negative means the model scored lower under dinucleotide matching than under GC matching.
    Almost every dataset is negative, which is the point: a harder negative set lowers the
    headline number even where it raises the measured contribution.
    """
    frame = cost.sort_values("cost")
    fig = go.Figure()
    fig.add_bar(
        x=frame["dataset"], y=frame["cost"],
        marker_color=[theme.BRASS if c < 0 else theme.STEEL for c in frame["cost"]],
        marker_line_width=0,
        customdata=frame[["protein", "cell"]].to_numpy(),
        hovertemplate="%{x}<br>%{customdata[0]} in %{customdata[1]}"
                      "<br>AUROC change %{y:+.4f}<extra></extra>",
        showlegend=False,
    )
    fig.add_hline(y=0, line_color=theme.INK_FAINT, line_width=1.2)
    lower = int((frame["cost"] < 0).sum())
    fig.add_annotation(
        xref="paper", yref="paper", x=0.02, y=0.06, showarrow=False, xanchor="left",
        text=f"{lower} of {len(frame)} datasets score LOWER under the stricter protocol",
        font={"size": 10.5, "color": theme.BRASS, "family": theme.MONO})
    fig.update_xaxes(showticklabels=False, title_text=f"{len(frame)} datasets, sorted")
    return _base(fig, 400, "Change in raw AUROC, dinucleotide minus GC")


def redraw_stability(draws: pd.DataFrame) -> go.Figure:
    """Five independent random draws of the same negative set. The result does not depend on one.

    An obvious objection is that the bias-aware number came from one lucky sample. It did not:
    five seeds land within 0.0009 of each other, against a protein-clustered standard error
    roughly four times that.
    """
    seeds = [c for c in draws["check"] if c.startswith("panel-mean contribution, seed")]
    vals = [draws[draws["check"] == s].iloc[0]["value"] for s in seeds]
    labels = [s.rsplit(" ", 1)[-1] for s in seeds]
    se = draws[draws["check"] == "between-protein standard error, the published draw"] \
        .iloc[0]["value"]
    mean = sum(vals) / len(vals)

    fig = go.Figure()
    fig.add_hrect(y0=mean - se, y1=mean + se, fillcolor=theme.STEEL, opacity=0.07,
                  line_width=0, layer="below")
    fig.add_hline(y=mean, line_color=theme.STEEL, line_width=1, line_dash="dot", opacity=0.7)
    fig.add_scatter(
        x=labels, y=vals, mode="markers+text",
        marker={"size": 13, "color": [theme.BRASS if i == 0 else theme.STEEL
                                      for i in range(len(vals))],
                "line": {"width": 2, "color": theme.SURFACE}},
        text=[f"{v:.5f}" for v in vals], textposition="top center",
        textfont={"size": 9.5, "color": theme.INK_FAINT},
        hovertemplate="seed %{x}<br>panel mean %{y:.5f}<extra></extra>", showlegend=False,
    )
    fig.add_annotation(
        xref="paper", yref="paper", x=0.99, y=0.06, showarrow=False, xanchor="right",
        text="shaded band = one protein-clustered standard error",
        font={"size": 10, "color": theme.INK_FAINT, "family": theme.MONO})
    span = max(vals) - min(vals)
    fig.update_yaxes(range=[mean - se * 1.7, mean + se * 1.7])
    fig.update_xaxes(title_text=f"random seed  ·  spread across draws {span:.5f}")
    return _base(fig, 400, "Panel-mean contribution, bias-aware arm")


def survey_timeline(per_method: pd.DataFrame) -> go.Figure:
    """Seven published methods, 2014 to 2025, and what each used for negatives.

    The chart is really about one column that is empty. Every surveyed method builds negatives
    somehow, and none of them reports what a composition-only baseline scores on the same data,
    so none of their headline numbers can be placed against a floor.

    Not a systematic review, and the view says so: seven sources chosen for coverage of the
    common constructions, not enumerated by a protocol.
    """
    kinds = {"coordinate": theme.STEEL, "other_rbp": theme.BRASS}
    frame = per_method.sort_values("year").copy()
    # Axis labels are shortened at a word boundary; the full name stays on hover. A 54-character
    # tick label reserves a quarter of the canvas and reads no better than a 24-character one.
    frame["short"] = [
        name if len(name) <= 26 else name[:26].rsplit(" ", 1)[0] + "..."
        for name in frame["method"]
    ]
    fig = go.Figure()
    for kind, colour in kinds.items():
        sub = frame[frame["kind"] == kind]
        if sub.empty:
            continue
        label = ("relocated genomic intervals" if kind == "coordinate"
                 else "other RBPs' binding sites")
        fig.add_scatter(
            x=sub["year"], y=sub["short"], mode="markers", name=label,
            marker={"size": 15, "color": colour, "symbol": "square",
                    "line": {"width": 2, "color": theme.SURFACE}},
            customdata=sub[["kind", "method"]].to_numpy(),
            hovertemplate="<b>%{customdata[1]}</b> (%{x})"
                          "<br>negatives: %{customdata[0]}<extra></extra>",
        )
    # Repeating "no composition baseline" against all seven points was noise, and the 2025
    # label ran off the axis. Said once, where it reads as the finding it is.
    fig.add_annotation(
        xref="paper", yref="paper", x=0.99, y=1.02, xanchor="right", yanchor="bottom",
        showarrow=False, text="NOT ONE REPORTS A COMPOSITION-ONLY BASELINE",
        font={"size": 10, "color": theme.BRASS, "family": theme.MONO})
    fig.update_xaxes(range=[2013, 2026.4], dtick=2)
    fig.update_yaxes(tickfont={"size": 10.5})
    return _base(fig, 420, "", "Year published")


def variant_ladder(ladder: pd.DataFrame) -> go.Figure:
    """The study this pivoted away from, as a ladder of controls.

    Read top down. Conservation alone, with no model at all, beats the model. Scoring each
    variant with a DIFFERENT protein's model still lands well above chance. The model is adding
    something real and far less than the headline suggests, and most of what looked
    protein-specific is not.
    """
    order = [
        ("conservation", "Conservation alone, no model", theme.SLATE),
        ("matched", "The protein's own model", theme.STEEL),
        ("mismatched", "A DIFFERENT protein's model", theme.BRASS),
        ("kmer", "4-mer composition only", theme.INK_FAINT),
    ]
    labels, vals, colours = [], [], []
    for arm, label, colour in order:
        row = ladder[ladder["arm"] == arm]
        if row.empty:
            continue
        labels.append(label)
        vals.append(float(row.iloc[0]["auroc"]))
        colours.append(colour)

    fig = go.Figure()
    fig.add_bar(
        y=labels, x=vals, orientation="h", marker_color=colours,
        marker_line_width=2, marker_line_color=theme.SURFACE,
        text=[f"{v:.3f}" for v in vals], textposition="outside",
        textfont={"size": 11.5, "color": theme.INK_MUTED},
        hovertemplate="%{y}<br>AUROC %{x:.4f}<extra></extra>", showlegend=False,
    )
    fig.add_vline(x=0.5, line_color=theme.INK_FAINT, line_width=1.4, line_dash="dot")
    fig.add_annotation(x=0.5, y=-0.62, text="chance", showarrow=False, yanchor="top",
                       font={"size": 10, "color": theme.INK_FAINT, "family": theme.MONO})
    _arrow(fig, vals[0], labels[0], "a control with no model at all<br>beats the model",
           ax=-96, ay=-46, colour=theme.SLATE)
    # Headroom on the right so the value labels and the callout both sit inside the frame.
    fig.update_xaxes(range=[0.45, 1.06])
    fig.update_yaxes(autorange="reversed")
    return _base(fig, 400, "", "AUROC, pathogenic vs benign non-coding variants")


def recommendation_intervals(rec: pd.DataFrame) -> go.Figure:
    """Whether the paper's own recommended rescaling survives correction. A forest plot.

    Drawn as intervals against zero because that is the question. The point estimate improves on
    all three protocol pairs, which is why the paper reports it; only one interval clears zero
    before correcting for the three comparisons, and none clears it after.
    """
    pairs = [("gc vs dn", "GC vs dinucleotide"), ("gc vs neg2", "GC vs bias-aware"),
             ("dn vs neg2", "Dinucleotide vs bias-aware")]
    fig = go.Figure()
    for style, suffix, colour, offset in (
        ("uncorrected", "", theme.STEEL, 0.16),
        ("Bonferroni, 3 pairs", ", Bonferroni over 3 pairs", theme.BRASS, -0.16),
    ):
        xs, ys, lo, hi = [], [], [], []
        for i, (key, _) in enumerate(pairs):
            row = rec[rec["check"] == f"rank agreement gain, {key}{suffix}"].iloc[0]
            xs.append(row["value"])
            ys.append(i + offset)
            lo.append(row["value"] - row["ci_low"])
            hi.append(row["ci_high"] - row["value"])
        fig.add_scatter(
            x=xs, y=ys, mode="markers", name=style,
            marker={"size": 11, "color": colour, "symbol": "diamond",
                    "line": {"width": 2, "color": theme.SURFACE}},
            error_x={"type": "data", "array": hi, "arrayminus": lo,
                     "color": colour, "thickness": 1.5, "width": 6},
            hovertemplate="%{x:+.4f}<extra>" + style + "</extra>",
        )
    # Zero is the decision line here, so it gets the emphasis a gridline would not carry.
    fig.add_vline(x=0, line_color=theme.INK, line_width=1.6, opacity=0.75)
    fig.add_annotation(
        x=0, y=2.62, text="NO IMPROVEMENT", showarrow=False, xanchor="center", yanchor="bottom",
        font={"size": 9, "color": theme.INK_MUTED, "family": theme.MONO})
    fig.update_yaxes(
        tickmode="array", tickvals=list(range(len(pairs))),
        ticktext=[label for _, label in pairs], range=[-0.6, 2.85],
        tickfont={"size": 11},
    )
    fig.update_xaxes(range=[-0.115, 0.185], tickformat="+.2f")
    return _base(fig, 400, "", "Gain in rank agreement from the recommended rescaling")


def class_ratio_robustness(ratios: pd.DataFrame) -> go.Figure:
    """The span at four negative-to-positive balances. A design choice, varied.

    One negative per positive is a convention, not a fact, and a result that only exists at 1:1
    would be an artefact of that convention. The span is drawn against 1.0, which is where it
    would sit if the protocol made no difference, and each point is annotated with whether the
    protocol ordering survived at that balance.
    """
    labels = ["1:1", "1:2", "1:4", "2:1"]
    order = ["2:1", "1:1", "1:2", "1:4"]        # negatives per positive, ascending
    spans, held, counts = [], [], []
    for label in order:
        spans.append(ratios[ratios["check"] == f"span, {label}"].iloc[0]["value"])
        held.append(ratios[ratios["check"] ==
                           f"ordering dn > gc > neg2 holds, {label}"].iloc[0]["value"])
        counts.append(int(ratios[ratios["check"] ==
                                 f"datasets in all three arms, {label}"].iloc[0]["value"]))
    assert set(order) == set(labels)

    fig = go.Figure()
    fig.add_scatter(
        x=order, y=spans, mode="lines+markers+text",
        line={"width": 2, "color": theme.STEEL, "shape": "spline"},
        marker={"size": 11, "color": [theme.STEEL if h == 1.0 else theme.BRASS for h in held],
                "line": {"width": 2, "color": theme.SURFACE}},
        text=[f"{s:.2f}x" for s in spans], textposition="top center",
        textfont={"size": 11, "color": theme.INK_MUTED},
        customdata=list(zip(counts, ["yes" if h == 1.0 else "NO" for h in held])),
        hovertemplate="%{x} negatives per positive<br>span %{y:.2f}x"
                      "<br>%{customdata[0]} datasets<br>ordering holds: %{customdata[1]}"
                      "<extra></extra>",
        showlegend=False,
    )
    fig.add_hline(y=1.0, line_color=theme.INK_MUTED, line_width=1.2, line_dash="dot")
    fig.add_annotation(
        xref="paper", x=0.01, y=1.0, text="no protocol effect", showarrow=False,
        xanchor="left", yanchor="bottom",
        font={"size": 10, "color": theme.INK_FAINT, "family": theme.MONO})
    published = ratios[ratios["check"] == "span at the published 1:1 balance"].iloc[0]["value"]
    _arrow(fig, "1:1", published, "the balance the paper reports", ax=0, ay=52,
           colour=theme.INK_FAINT)
    fig.update_yaxes(range=[0, max(spans) * 1.28])
    return _base(fig, 400, "Span across the three protocols",
                 "Negatives per positive")


def protocol_trajectories(per_dataset: pd.DataFrame) -> go.Figure:
    """Every dataset as one line across the three protocols. The heatmap, made readable.

    A 94-row heatmap gives each dataset a one-pixel stripe and most of the values sit near zero
    against a maximum of 0.24, so it renders as a dark wall with no structure visible. The same
    data as trajectories shows what the heatmap was meant to: the lines move together, which
    means the protocol sets the level and the dataset only shifts it.

    Lines are coloured by whether that dataset follows the panel ordering. The panel median is
    drawn over the top, thick, because 94 thin lines need an anchor.
    """
    arms = ("gc", "dn", "neg2")
    xs = [theme.PROTOCOL_LABEL[a] for a in arms]
    frame = per_dataset.dropna(subset=[f"gain_{a}" for a in arms])

    follows = 0
    fig = go.Figure()
    for _, row in frame.iterrows():
        ys = [row[f"gain_{a}"] for a in arms]
        # The panel ordering is dinucleotide highest, bias-aware lowest.
        agrees = ys[1] > ys[0] > ys[2]
        follows += agrees
        fig.add_scatter(
            x=xs, y=ys, mode="lines",
            line={"width": 1, "color": theme.STEEL if agrees else theme.BRASS,
                  "shape": "spline", "smoothing": 0.6},
            opacity=0.34 if agrees else 0.5,
            customdata=[[row["dataset"]]] * 3,
            hovertemplate="%{customdata[0]}<br>%{x}  %{y:.4f}<extra></extra>",
            showlegend=False,
        )

    medians = [frame[f"gain_{a}"].median() for a in arms]
    fig.add_scatter(
        x=xs, y=medians, mode="lines+markers", name="panel median",
        line={"width": 3, "color": theme.INK, "shape": "spline", "smoothing": 0.6},
        marker={"size": 11, "color": theme.INK,
                "line": {"width": 3, "color": theme.SURFACE}},
        hovertemplate="panel median<br>%{x}  %{y:.4f}<extra></extra>", showlegend=False,
    )
    for x, value in zip(xs, medians):
        fig.add_annotation(
            x=x, y=value, text=f"{value:.4f}", showarrow=False, yshift=20,
            font={"size": 11, "color": theme.INK, "family": theme.MONO},
            bgcolor="rgba(16,20,28,0.85)", borderpad=3)

    fig.add_hline(y=0, line_color=theme.BORDER, line_width=1)
    fig.add_annotation(
        xref="paper", yref="paper", x=0.01, y=0.99, showarrow=False, xanchor="left",
        align="left",
        text=(f"<span style='color:{theme.STEEL}'>{follows} of {len(frame)} datasets</span>"
              " follow the panel ordering<br>"
              f"<span style='color:{theme.BRASS}'>{len(frame) - follows}</span> do not"),
        font={"size": 11, "family": theme.MONO, "color": theme.INK_MUTED})
    return _base(fig, 520, "Nested contribution (AUROC points)")
