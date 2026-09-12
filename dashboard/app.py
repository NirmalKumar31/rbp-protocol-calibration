"""RBP Benchmark Calibration Explorer.

A reader over the committed tables in `results/tables/`. It opens CSV files, filters them in
memory and draws them. It does not run an analysis, fit a model, or write to the repository, and
it makes no network request for data. The browser does fetch two font families from Google
Fonts on first paint; see "Fonts and network" in dashboard/README.md.

One rule shapes the whole interface. A number that came from a committed table is marked
PUBLISHED and carries its bootstrap interval. A number computed live from whatever the user
filtered to is marked LIVE, carries no interval, and is never called an estimate. Filtering is
the interactive part; the published estimand does not move when you filter.

Run with `streamlit run dashboard/app.py` from the repository root.
"""

from __future__ import annotations

import streamlit as st
from rbp_dashboard import copy, data, figures, graphics, theme

_NO_SHARED_SPAN = "SPAN excluding every protein shared with our panel"

PAPER_DOI = "https://doi.org/10.20944/preprints202609.0883.v1"
CODE_DOI = "https://doi.org/10.5281/zenodo.22679284"

st.set_page_config(
    page_title="RBP Benchmark Calibration Explorer",
    page_icon=None,
    layout="wide",
    # "auto", not "expanded". Streamlit collapses an auto sidebar below its mobile breakpoint;
    # forcing it open put a 330px navigation panel over a 430px viewport and clipped the body
    # text mid-sentence. On desktop "auto" still opens it.
    initial_sidebar_state="auto",
)
st.markdown(theme.CSS, unsafe_allow_html=True)
st.markdown(theme.background(), unsafe_allow_html=True)

CHART = {"displayModeBar": False, "scrollZoom": False}
CHART_ZOOM = {"displaylogo": False, "scrollZoom": True,
              "modeBarButtonsToRemove": ["lasso2d", "select2d", "autoScale2d"]}


# ------------------------------------------------------------------------------------ helpers

def src(*names: str) -> None:
    st.markdown(
        '<p class="src">' + "  ·  ".join(f"results/tables/{n}" for n in names) + "</p>",
        unsafe_allow_html=True,
    )


def keys(*items: tuple[str, str]) -> None:
    """A band of the figures that carry the argument, set larger than body text."""
    cells = "".join(
        f'<div><div class="key-v">{value}</div><div class="key-l">{label}</div></div>'
        for value, label in items
    )
    st.markdown(f'<div class="keys">{cells}</div>', unsafe_allow_html=True)


def pull(text: str) -> None:
    st.markdown(f'<div class="pull">{text}</div>', unsafe_allow_html=True)


def diagram(svg: str, caption: str) -> None:
    """A schematic, drawn as page content. Not data, and labelled as such inside the SVG."""
    st.markdown(f'<div class="fig">{svg}</div>', unsafe_allow_html=True)
    st.markdown(f'<p class="plain">{caption}</p>', unsafe_allow_html=True)


def chart(fig, key: str, *sources: str, zoom: bool = False) -> None:
    """Draw a figure with its plain-language line and its source files. All three or none.

    Routed through one function so a chart cannot reach the page without a sentence explaining
    it to someone outside the field, and without naming the file it came from. Both are
    positional arguments rather than optional keywords for that reason.
    """
    # theme=None, not the default "streamlit". Streamlit's own theme REPLACES the figure
    # template rather than merging with it, which rendered every chart on a white paper while
    # the template asked for a transparent one.
    st.plotly_chart(fig, use_container_width=True, theme=None,
                    config=CHART_ZOOM if zoom else CHART)
    if key not in copy.PLAIN:
        raise KeyError(f"no plain-language caption for {key!r}")
    st.markdown(f'<p class="plain">{copy.PLAIN[key]}</p>', unsafe_allow_html=True)
    if sources:
        src(*sources)


def box(kind: str, text: str) -> None:
    st.markdown(f'<div class="{kind}">{text}</div>', unsafe_allow_html=True)


def badge(kind: str) -> str:
    return ('<span class="pub">published</span>' if kind == "pub"
            else '<span class="live">live · descriptive</span>')


# Reset by page(). Section numbers are the reading order, which is otherwise invisible on a
# page this long, and numbering them by hand at thirty call sites would drift the first time a
# section moved.
_SECTION = {"n": 0}


def heading(text: str, kind: str | None = None) -> None:
    _SECTION["n"] += 1
    st.markdown(
        f'<div class="step"><span class="step-n">{_SECTION["n"]:02d}</span>'
        f'<span class="step-t">{text}</span>'
        f'{badge(kind) if kind else ""}</div>',
        unsafe_allow_html=True,
    )


def page(title: str, guide: str | None = None) -> None:
    """Title, hairline rule, and the reading guide. Every view opens the same way.

    The guide exists because the first build was legible to someone who already knew the result
    and opaque to everyone else. `filters` is stated on every view because four of seven show
    only panel-level published estimates, which correctly do not move when a protein is
    deselected, and silence about that reads as a broken filter.
    """
    _SECTION["n"] = 0
    st.markdown(f"# {title}")
    st.markdown(theme.RULE, unsafe_allow_html=True)
    g = copy.GUIDE.get(guide or title)
    if not g:
        return
    label, explain = copy.FILTER_STATE[g["filters"]]
    chip = f'<span class="live">{label}</span>' if label else ""
    tail = (f'<p style="margin:0.8rem 0 0 0;font-size:0.79rem;color:#616D7E">{explain}</p>'
            if explain else "")
    items = "".join(f"<li>{c}</li>" for c in g["checks"])
    st.markdown(
        '<div class="guide">'
        f'<div class="gl">What this shows{chip}</div>'
        f"<p>{g['lead']}</p>"
        f'<div class="gl">What to check, in order</div><ol>{items}</ol>'
        f"{tail}</div>",
        unsafe_allow_html=True,
    )
    if (guide or title) in copy.NARRATOR:
        with st.expander("Read the background", expanded=False):
            st.markdown(copy.NARRATOR[guide or title])


def term(name: str) -> str:
    return copy.GLOSSARY[name]


def guard(*tables: str) -> bool:
    try:
        for name in tables:
            data.load(name)
    except data.SchemaError as exc:
        st.error(f"This view cannot render: {exc}")
        return False
    return True


NO_FILTER = {"cells": None, "proteins": None, "size_range": None, "size_column": None}


def sidebar_filters(frame, size_column: str, active: bool):
    """Global filters, rendered only on views that use them.

    Offering a control that does nothing is worse than offering none: the reader filters, the
    page does not move, and the reasonable conclusion is that the dashboard is broken. Seven of
    the ten views never call `narrow()`, so on those the controls are not drawn at all and this
    returns a filter that matches everything.
    """
    if not active:
        st.markdown(
            '<div class="sidebar-title">Filter</div>'
            '<p class="src" style="margin:0">not used on this view.<br>'
            "open Dataset explorer, Protocol sensitivity<br>or Model comparison to filter.</p>",
            unsafe_allow_html=True,
        )
        return NO_FILTER
    facets = data.facets(frame)
    st.markdown('<div class="sidebar-title">Filter</div>', unsafe_allow_html=True)
    st.markdown(
        '<p class="src" style="margin:0 0 0.7rem 0">narrows per-dataset charts only.<br>'
        "published panel estimates are fixed.</p>",
        unsafe_allow_html=True,
    )
    cells = st.multiselect("Cell line", facets["cells"], default=facets["cells"])
    proteins = st.multiselect(
        "Protein", facets["proteins"], default=[],
        help="Empty means all proteins. Filtering here narrows every per-dataset chart.",
    )
    lo, hi = int(frame[size_column].min()), int(frame[size_column].max())
    size = st.slider(
        "Windows per dataset", lo, hi, (lo, hi), step=max(1, (hi - lo) // 100),
        help="The number of scored windows. Small datasets carry wider dataset-level noise.",
    )
    return {"cells": cells, "proteins": proteins, "size_range": size,
            "size_column": size_column}


def narrow(frame, flt):
    return data.apply_filters(
        frame, cells=flt["cells"], proteins=flt["proteins"],
        size_range=flt["size_range"], size_column=flt["size_column"],
    )


def filter_readout(full, subset) -> None:
    if len(subset) == len(full):
        st.markdown(
            f'<p class="src">all {len(full)} datasets · no filter active</p>',
            unsafe_allow_html=True)
    else:
        st.markdown(
            f'<p class="src">{len(subset)} of {len(full)} datasets after filtering '
            f'· published estimates below are unchanged and still cover all {len(full)}</p>',
            unsafe_allow_html=True)


# -------------------------------------------------------------------------------------- views

def view_overview(flt) -> None:
    page("RBP Benchmark Calibration Explorer", guide="Overview")
    st.markdown(
        "How much does an RNA sequence model add beyond plain nucleotide composition, and how "
        "much does that answer depend on how the *unbound* windows were chosen?"
    )
    diagram(
        graphics.hero(),
        "A protein sits on one stretch of an RNA strand. That stretch is a <b>bound window</b>. "
        "Every other stretch is a candidate for the <b>unbound</b> comparison set, and which "
        "ones you choose is what this study measures.",
    )
    if not guard("cross_fitting.csv", "three_arm_contrast.csv", "panel_summary.csv"):
        return

    panel = data.table("panel_summary.csv")
    span, lo, hi = data.interval("cross_fitting.csv", "4-mer three-arm span, fully cross-fitted")
    two_stage = data.value("cross_fitting.csv", "4-mer three-arm span, as published")
    removed = min(
        data.value("cross_fitting.csv",
                   f"2-mer fraction of the published value removed by cross-fitting, {a} arm")
        for a in ("gc", "dn", "neg2")
    )

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Primary span", f"{span:.2f}x",
              help="Largest over smallest nested contribution across the three negative-set "
                   "protocols, 4-mer logistic regression, cross-fitted estimator. " + term("span"))
    # No `or 94` fallback here. It was hiding the case where the column stops parsing as the
    # string "True", which would silently print the right number from the wrong source.
    if panel["in_both_arms"].dtype == bool:
        in_all_three = int(panel["in_both_arms"].sum())
    else:
        in_all_three = int((panel["in_both_arms"].astype(str) == "True").sum())
    if in_all_three == 0:
        raise data.SchemaError(
            "panel_summary.csv: in_both_arms parsed to zero datasets, so its values are not "
            "the booleans this count assumes")
    c2.metric("Datasets", f"{in_all_three}", help=copy.PANEL_CAVEAT)
    c3.metric("Protocols", "3", help="GC-matched, dinucleotide-matched, bias-aware.")
    c4.metric("Known-null bias removed", f">={removed * 100:.0f}%",
              help="On a constructed test whose true answer is zero, cross-fitting removes at "
                   "least this share of what the conventional estimator reports.")
    st.markdown(
        f'<p class="src">95% CI {lo:.2f} to {hi:.2f}, protein-clustered bootstrap  ·  '
        f'two-stage estimator gives {two_stage:.2f}x on the same data</p>',
        unsafe_allow_html=True)

    pull(
        "Three ways of choosing the comparison sequences, all of them in the literature, with "
        "the model class, source peaks, fold design and estimator held fixed. The measured "
        f"answer moves <b>{span:.2f}-fold</b> between them."
    )

    heading("The finding", "pub")
    st.markdown(
        "Hold the model class, the source peaks, the chromosome-blocked folds and the estimator "
        f"fixed. Change how the negative windows are built. The contribution moves **{span:.2f}"
        "-fold**, and it moves *against* the apparent score: the protocol that looks easiest by "
        "headline AUROC is where the model contributes least."
    )
    chart(
        figures.protocol_contribution(data.table("three_arm_contrast.csv")),
        "protocol_contribution", "three_arm_contrast.csv")

    box("note", copy.NOT_ONLY_NEGATIVES)
    box("warn", "<strong>What this is.</strong> " + copy.NOT_MONITORING)

    heading("Where to look")
    a, b = st.columns(2)
    with a:
        st.markdown(
            "- **Protocol sensitivity** — the three protocols, paired per dataset\n"
            "- **Model comparison** — does it survive a change of model class\n"
            "- **Cross-fitting** — the estimator tested against a known zero"
        )
    with b:
        st.markdown(
            "- **Dataset explorer** — filter, rank and drill into single datasets\n"
            "- **External validation** — what replicated, and what did not\n"
            "- **Reproducibility** — provenance of every committed table"
        )
    st.markdown(
        f'<p class="src">preprint <a href="{PAPER_DOI}">10.20944/preprints202609.0883.v1</a>'
        f"  ·  code and evidence <a href='{CODE_DOI}'>10.5281/zenodo.22679284</a></p>",
        unsafe_allow_html=True)




def view_why(flt) -> None:
    page("Why this exists")
    if not guard("negative_set_survey.csv", "negative_set_survey_per_method.csv",
                 "variant_ladder.csv"):
        return

    survey = data.table("negative_set_survey.csv")
    n_methods = int(data.value("negative_set_survey.csv", "methods and benchmarks surveyed"))
    n_baseline = int(data.value("negative_set_survey.csv",
                                "surveyed sources reporting a composition-only baseline"))

    keys(
        (f"{n_baseline}<span class='key-s'> of {n_methods}</span>",
         "report a composition baseline"),
        ("2014&ndash;2025", "years covered"),
        ("7", "constructions surveyed"),
    )

    box("warn",
        f"<strong>{n_baseline} of {n_methods}.</strong> Every surveyed method builds a negative "
        "set, and not one reports what a composition-only model scores on that same negative "
        "set. Without that floor a headline AUROC cannot be placed against anything.")

    heading("What the field reports", "pub")
    chart(figures.survey_timeline(data.table("negative_set_survey_per_method.csv")),
          "survey_timeline", "negative_set_survey_per_method.csv")
    box("note",
        "This is a <strong>targeted survey, not a systematic review</strong>. Seven sources were "
        "chosen to cover the common negative-set constructions, not enumerated under a search "
        "protocol, and the paper describes it that way.")

    with st.expander("The seven sources, with quotes"):
        per = data.table("negative_set_survey_per_method.csv")
        for _, row in per.sort_values("year").iterrows():
            st.markdown(
                f"**{row['method']}** ({int(row['year'])}) · negatives: `{row['kind']}`  \n"
                f"> {row['quote']}  \n"
                f"[source]({row['url']})"
            )
        src("negative_set_survey_per_method.csv")

    heading("Where this came from", "pub")
    st.markdown(
        "Before this, I was asking whether these models could flag disease-causing mutations. "
        "This is the ladder of controls that ended that line of work."
    )
    chart(figures.variant_ladder(data.table("variant_ladder.csv")),
          "variant_ladder", "variant_ladder.csv")
    box("bad",
        "<strong>Conservation alone scores 0.908; the model scores 0.829.</strong> A control "
        "with no model in it beats the model. And a <em>different</em> protein's model still "
        "reaches 0.680 against a chance level of 0.500, so most of what looked protein-specific "
        "was not.")
    st.markdown(
        "The model's contribution survived adjustment, but it was much smaller than the "
        "headline suggested. Chasing a better number was the wrong response. Measuring what the "
        "number rests on was the right one, and that is the rest of this site."
    )

    with st.expander("Survey counts"):
        st.dataframe(survey, use_container_width=True, hide_index=True)
        src("negative_set_survey.csv")


def view_recommendation(flt) -> None:
    page("Does the fix work?")
    if not guard("recommendation_works.csv"):
        return

    rec = data.table("recommendation_works.csv")
    improves = int(data.value("recommendation_works.csv",
                              "protocol pairs where rank agreement improves"))
    total = int(data.value("recommendation_works.csv", "number of protocol pairs"))
    shrinks = int(data.value("recommendation_works.csv",
                             "protocol pairs where disagreement shrinks"))

    c1, c2, c3 = st.columns(3)
    c1.metric("Pairs where it improves", f"{improves} of {total}")
    c2.metric("Pairs where disagreement shrinks", f"{shrinks} of {total}")
    c3.metric("Pairs clearing zero after correction", f"0 of {total}")

    heading("The pre-specified test", "pub")
    st.markdown(
        "If rescaling by headroom works, two datasets ranked differently under two protocols "
        "should agree more once rescaled. That is measurable, and I wrote down how I would "
        "measure it before looking."
    )
    chart(figures.recommendation_intervals(rec), "recommendation_intervals",
          "recommendation_works.csv")

    box("bad",
        "<strong>The recommendation does not survive its own test.</strong> All three point "
        "estimates are positive, so the direction is consistent. But only one of three intervals "
        "clears zero before correcting for three comparisons, and after Bonferroni correction "
        "none of them does.")
    gain, raw_lo, raw_hi = data.interval("recommendation_works.csv",
                                        "rank agreement gain, gc vs dn")
    _, bon_lo, bon_hi = data.interval(
        "recommendation_works.csv", "rank agreement gain, gc vs dn, Bonferroni over 3 pairs")
    st.markdown(
        f"The strongest pair, GC against dinucleotide, gains **{gain:+.4f}**. Uncorrected its "
        f"interval is [{raw_lo:+.4f}, {raw_hi:+.4f}], which excludes zero. Corrected for the "
        f"three comparisons it is [{bon_lo:+.4f}, {bon_hi:+.4f}], which does not. The paper "
        "reports the corrected one."
    )
    box("note",
        "This view exists because leaving it out would make the paper tidier and less true. A "
        "recommendation that fails its own pre-specified test is still a finding, and the "
        "alternative is quoting the uncorrected interval and hoping nobody checks.")

    with st.expander("All twenty-one measured quantities"):
        st.dataframe(rec, use_container_width=True, hide_index=True)
        src("recommendation_works.csv")


def view_construction(flt) -> None:
    page("How negatives are built")
    if not guard("match_quality.csv", "cost_of_matching.csv", "negative_draws.csv",
                 "class_ratio.csv"):
        return

    match = data.table("match_quality.csv")

    diagram(
        graphics.matching_procedure(),
        "Take a window the protein binds, measure it, and look for an unbound window to pair with "
        "it. The two composition-matched protocols search free genomic intervals and differ in "
        "how strict the match must be. <b>The bias-aware protocol is not a third strictness "
        "setting</b>: it draws from a different pool entirely, other proteins' binding sites, "
        "and matches fold rather than composition. The bars are how often each ends up "
        "composition-matched anyway.",
    )

    heading("How closely each protocol actually matched", "pub")
    st.markdown(
        "Every bound window is paired with an unbound one. This is how similar those pairs "
        "turned out to be, measured on the committed evidence rather than asserted."
    )
    chart(figures.animated_matching(match), "animated_matching", "match_quality.csv")

    with st.expander("All three curves at once"):
        chart(figures.match_quality_curve(match), "match_quality_curve", "match_quality.csv")

    c1, c2, c3 = st.columns(3)
    for col, arm in zip((c1, c2, c3), ("gc", "dn", "neg2")):
        frac = data.value("match_quality.csv",
                          f"fraction of pairs within |dGC| 0.05, {arm} arm")
        med = data.value("match_quality.csv", f"gc_gap_median, {arm} arm")
        col.metric(theme.PROTOCOL_LABEL[arm], f"{frac:.0%}",
                   help=f"Share of pairs whose GC fraction differs by less than 0.05. "
                        f"Median gap {med:.4f}.")
    st.markdown(
        '<p class="src">share of bound/unbound pairs matched within 0.05 GC fraction</p>',
        unsafe_allow_html=True)

    box("note",
        "This is the mechanism the rest of the study measures. The bias-aware protocol matches "
        "only <strong>20%</strong> of its pairs on composition, against 95% for GC-matched, "
        "because it is not trying to. It corrects for the assay's technical bias instead. That "
        "leaves a composition baseline that already scores 0.825, and a baseline that high "
        "leaves little room for a model to add anything on top.")

    heading("The gap itself", "pub")
    chart(figures.match_gap_bars(match), "match_gap_bars", "match_quality.csv")
    st.markdown(
        "The dinucleotide protocol achieves a *tighter* median GC gap than the GC protocol "
        f"({data.value('match_quality.csv', 'gc_gap_median, dn arm'):.4f} against "
        f"{data.value('match_quality.csv', 'gc_gap_median, gc arm'):.4f}) while also cutting "
        "the dinucleotide mismatch by "
        f"{data.value('match_quality.csv', 'dinucleotide L1 improvement factor, dn vs gc arm'):.2f}"
        "-fold. It is strictly the harder constraint."
    )

    heading("What stricter matching costs", "pub")
    st.markdown(
        "Harder negatives lower the model's raw score. That is expected. The finding is that "
        "the model's *measured contribution* moves the other way."
    )
    chart(figures.matching_cost(data.table("cost_of_matching.csv")),
          "matching_cost", "cost_of_matching.csv")

    heading("Does it survive the design choices?", "pub")
    st.markdown(
        "Pairing one unbound window with each bound one is a convention, not a fact. If the "
        "result only existed at that balance it would be an artefact of the convention."
    )
    chart(figures.class_ratio_robustness(data.table("class_ratio.csv")),
          "class_ratio_robustness", "class_ratio.csv")
    box("good",
        "The span ranges from "
        f"{data.value('class_ratio.csv', 'span, 1:1'):.2f}x to "
        f"{data.value('class_ratio.csv', 'span, 1:4'):.2f}x across the four balances tested, "
        "and the protocol ordering holds at every one of them. Changing the balance moves the "
        "number a little and removes nothing.")

    heading("Is it just one lucky draw?", "pub")
    chart(figures.redraw_stability(data.table("negative_draws.csv")),
          "redraw_stability", "negative_draws.csv")
    between_protein = data.value("negative_draws.csv",
                                 "between-protein standard error, the published draw")
    width_ratio = data.value("negative_draws.csv",
                             "ratio of combined to published interval width")
    st.markdown(
        "Five independent draws of the bias-aware negative set. The spread between them is "
        f"{data.value('negative_draws.csv', 'between-draw SD of the panel mean'):.5f}, against "
        "a protein-clustered standard error of "
        f"{between_protein:.5f}. "
        "Folding draw uncertainty in would widen the published interval by "
        f"{(width_ratio - 1) * 100:.1f}%."
    )

    with st.expander("Match quality, full table"):
        st.dataframe(match, use_container_width=True, hide_index=True)
        src("match_quality.csv")


def view_protocol(flt) -> None:
    page("Protocol sensitivity")
    if not guard("three_arm_contrast.csv", "three_arm_per_dataset.csv"):
        return

    contrast = data.table("three_arm_contrast.csv")
    full = data.table("three_arm_per_dataset.csv")
    subset = narrow(full, flt)
    filter_readout(full, subset)

    diagram(
        graphics.negative_sets(),
        "The top line is the transcript with the windows the protein actually binds. Each row "
        "below is one procedure for choosing the comparison windows. The first two search free "
        "genomic intervals at different strictness; the third draws from other proteins' "
        "binding sites instead. The measured answer moves almost fivefold between them.",
    )

    with st.expander("What the three protocols are"):
        for key in ("GC-matched", "dinucleotide-matched", "bias-aware"):
            st.markdown(f"**{key}.** {term(key)}")

    tab1, tab2, tab3 = st.tabs(
        ["1 · Levels", "2 · Per dataset", "3 · The inverse relation"])

    with tab1:
        diagram(
            graphics.composition(),
            "The thing the model has to beat. It counts how often each letter and each adjacent "
            "pair of letters appears, and nothing else: no order, no motifs, no learning.",
        )

        heading("Baseline and apparent AUROC", "pub")
        st.markdown(
            "Two quantities on one scale. The first is what plain composition scores alone; the "
            "second is composition plus the model's score, which is what a paper prints."
        )
        chart(
            figures.protocol_levels(contrast),
            "protocol_levels", "three_arm_contrast.csv")

        heading("What the model adds", "pub")
        chart(
            figures.protocol_contribution(contrast),
            "protocol_contribution", "three_arm_contrast.csv")
        box("note",
            "Read them together. The bias-aware protocol has the <em>highest</em> apparent AUROC "
            "and the <em>lowest</em> contribution. Dinucleotide-matched is the reverse. A headline "
            "AUROC alone cannot distinguish them.")

    with tab2:
        if subset.empty:
            st.warning("No datasets match the current filter.")
        else:
            heading("Distribution", "live")
            chart(
                figures.animated_protocol_distribution(subset),
                "animated_protocol_distribution", "three_arm_per_dataset.csv")
            with st.expander("All three at once, as boxes"):
                chart(
                    figures.contribution_distribution(subset),
                    "contribution_distribution", "three_arm_per_dataset.csv", zoom=True)

            heading("Paired within dataset", "live")
            left, _ = st.columns([1, 2])
            with left:
                pair = st.selectbox("Compare", ["dn vs gc", "neg2 vs gc", "dn vs neg2"])
            mapping = {"dn vs gc": ("gc", "dn"), "neg2 vs gc": ("gc", "neg2"),
                       "dn vs neg2": ("neg2", "dn")}
            a, b = mapping[pair]
            chart(
                figures.paired_slopes(subset, a, b, limit=len(subset)),
                "paired_slopes", "three_arm_per_dataset.csv", zoom=True)

            heading("Filtered subset against the published estimate")
            arm = st.radio("Protocol", ["dn", "gc", "neg2"], horizontal=True,
                           format_func=lambda a: theme.PROTOCOL_LABEL[a])
            stat = data.descriptive(subset, f"gain_{arm}")
            published = data.interval("three_arm_contrast.csv",
                                      f"nested contribution, {arm} arm")
            chart(
                figures.filtered_versus_published(stat, published, arm),
                "filtered_versus_published",
                "three_arm_contrast.csv", "three_arm_per_dataset.csv")
            box("warn",
                "The hatched bar is the mean of the datasets you filtered to. It is "
                "<strong>descriptive only</strong>: no protein-clustered interval, no estimand, "
                "no claim. The solid bar is the published panel estimate over all 94 datasets "
                "and does not move when you filter.")

    with tab3:
        heading("Apparent AUROC against contribution", "live")
        st.markdown(
            "Within this panel the two move in opposite directions. **This is the finding that "
            "did not replicate on the held-out benchmark.**"
        )
        if subset.empty:
            st.warning("No datasets match the current filter.")
        else:
            chart(
                figures.apparent_versus_contribution(subset),
                "apparent_versus_contribution", "three_arm_per_dataset.csv", zoom=True)

    with st.expander("Numbers behind these charts"):
        st.dataframe(contrast[["check", "value", "ci_low", "ci_high", "n"]],
                     use_container_width=True, hide_index=True)
        src("three_arm_contrast.csv")


def view_models(flt) -> None:
    page("Model comparison")
    if not guard("three_arm_models.csv", "three_arm_models_per_dataset.csv",
                 "cross_fitting.csv"):
        return

    models = data.table("three_arm_models.csv")
    box("warn", "<strong>Estimator notice.</strong> " + copy.ESTIMATOR_CAVEAT)

    heading("Protocol ordering, model by model", "pub")
    st.markdown(
        "Press **PLAY**, or step with the slider. The bars move; the ordering does not. All three "
        "model classes put their lowest contribution in the bias-aware arm."
    )
    chart(
        figures.animated_model_walk(models),
        "animated_model_walk", "three_arm_models.csv")

    heading("How far each model's answer moves", "pub")
    chart(
        figures.model_spans(models),
        "model_spans", "three_arm_models.csv")
    st.markdown(
        "All three spans sit far from 1.0, so protocol dependence is not an artefact of one model "
        "class. They are **not** comparable as estimates of the same quantity: the 4-mer figure "
        "here is two-stage so the three sit on one footing, while the paper's primary 4-mer span "
        f"is the cross-fitted "
        f"{data.value('cross_fitting.csv', '4-mer three-arm span, fully cross-fitted'):.2f}x."
    )

    heading("Side by side", "pub")
    chart(
        figures.model_by_protocol(models),
        "model_by_protocol", "three_arm_models.csv")

    with st.expander("Per-dataset values"):
        st.dataframe(narrow(data.table("three_arm_models_per_dataset.csv"), flt),
                     use_container_width=True, hide_index=True)
        src("three_arm_models_per_dataset.csv")


def view_explorer(flt) -> None:
    page("Dataset explorer")
    if not guard("three_arm_per_dataset.csv", "three_arm_models_per_dataset.csv",
                 "panel_summary.csv"):
        return

    full = data.table("three_arm_per_dataset.csv")
    subset = narrow(full, flt)
    filter_readout(full, subset)
    if subset.empty:
        st.warning("No datasets match the current filter. Widen it in the sidebar.")
        return

    tab1, tab2, tab3 = st.tabs(
        ["1 · Map", "2 · Ranking", "3 · Single dataset"])

    with tab1:
        heading("Every dataset across every protocol", "live")
        st.markdown("One line per dataset. Hover any line for its name and its three values.")
        chart(
            figures.protocol_trajectories(subset),
            "protocol_trajectories", "three_arm_per_dataset.csv", zoom=True)

        with st.expander("The same data as a heatmap"):
            sort_by = st.radio("Sort by", ["dn", "gc", "neg2"], horizontal=True,
                               format_func=lambda a: theme.PROTOCOL_LABEL[a])
            chart(
                figures.dataset_heatmap(subset, sort_by),
                "dataset_heatmap", "three_arm_per_dataset.csv", zoom=True)
        box("note",
            "Read this carefully, because it is easy to overstate. The shift is systematic at "
            "the level of the panel mean, which is what moves 4.84-fold. It is "
            "<strong>not</strong> uniform across datasets: 59 of 94 follow the panel ordering "
            "and 35 do not. Nor does the protocol explain most of the variation, because "
            "dataset-to-dataset spread is the larger of the two, at about 1.6 times the "
            "between-protocol spread. Both hold at once: the claim is about where the panel "
            "mean lands, not about which source of variance is bigger.")

    with tab2:
        heading("Strongest datasets", "live")
        c1, c2 = st.columns([1, 1])
        with c1:
            arm = st.radio("Protocol", ["dn", "gc", "neg2"], horizontal=True, key="rank_arm",
                           format_func=lambda a: theme.PROTOCOL_LABEL[a])
        with c2:
            # A slider needs a range. Filter to five or fewer datasets and min_value would meet
            # or exceed max_value, which Streamlit raises on rather than clamping. This crashed
            # in the browser at one dataset; the first harness only ever tested 49 and 94.
            if len(subset) > 6:
                top = st.slider("How many to show", 5, min(40, len(subset)),
                                min(20, len(subset)))
            else:
                top = len(subset)
                st.markdown(
                    f'<p class="src" style="margin-top:1.9rem">showing all {top} '
                    "datasets in the filter</p>", unsafe_allow_html=True)
        chart(
            figures.ranked_datasets(subset, arm, top),
            "ranked_datasets", "three_arm_per_dataset.csv")

        stat = data.descriptive(subset, f"gain_{arm}")
        if stat:
            m1, m2, m3, m4 = st.columns(4)
            m1.metric("Datasets", stat["n"])
            m2.metric("Mean", f"{stat['mean']:.4f}")
            m3.metric("Median", f"{stat['median']:.4f}")
            m4.metric("Range", f"{stat['min']:.3f} to {stat['max']:.3f}")
            st.markdown(
                f'<p class="src">descriptive statistics of the filtered subset, no interval  ·  '
                f'the published panel estimate for '
                f'{theme.PROTOCOL_LABEL[arm]} is '
                f'{data.ci_text("three_arm_contrast.csv", f"nested contribution, {arm} arm", 4)}'
                "</p>", unsafe_allow_html=True)

    with tab3:
        heading("Drill into one dataset", "live")
        models_per = data.table("three_arm_models_per_dataset.csv")
        options = sorted(set(subset["dataset"]) & set(models_per["dataset"]))
        if not options:
            st.warning("No dataset in the filter has model-level values.")
            return
        chosen = st.selectbox("Dataset", options,
                              help="Filtered by the sidebar. Type to search.")
        row3 = subset[subset["dataset"] == chosen].iloc[0]
        rowm = models_per[models_per["dataset"] == chosen].iloc[0]

        m1, m2, m3 = st.columns(3)
        m1.metric("Protein", str(row3["protein"]))
        m2.metric("Cell line", str(row3["cell"]))
        m3.metric("Windows", f"{int(row3['n_dn']):,}")

        chart(
            figures.dataset_profile(rowm),
            "dataset_profile", "three_arm_models_per_dataset.csv")
        box("warn",
            "Single-dataset values, shown because they are in the committed table. They are "
            "<strong>not</strong> estimates: the study's unit of inference is the panel with a "
            "protein-clustered interval, and one dataset carries no interval at all. Small "
            "datasets in particular move a great deal on resampling.")

        with st.expander("This dataset's three-arm row"):
            st.dataframe(subset[subset["dataset"] == chosen], use_container_width=True,
                         hide_index=True)
            src("three_arm_per_dataset.csv")


def view_crossfitting(flt) -> None:
    page("Cross-fitting")
    st.markdown("What happens when the estimator is tested on a question whose answer is known.")
    if not guard("cross_fitting.csv", "estimator_floor.csv"):
        return

    cross = data.table("cross_fitting.csv")

    with st.expander("The two estimators, in plain language", expanded=True):
        st.markdown(f"**Two-stage.** {term('two-stage estimator')}")
        st.markdown(f"**Cross-fitted.** {term('cross-fitted estimator')}")
        st.markdown(f"**Estimator floor.** {term('estimator floor')}")

    diagram(
        graphics.cross_fitting(),
        "Five slices of the same data. In the top route the score being judged was produced by "
        "a model that had already seen the held-out slice, so information leaks backwards. In "
        "the bottom route it had not. That one difference is what the next chart measures.",
    )

    heading("The known-answer test", "pub")
    st.markdown(
        "Score each window with a **2-mer** model. The composition baseline already contains "
        "every 2-mer frequency, so a 2-mer model can add nothing: the true contribution is "
        "**zero by construction**. Anything reported here is the estimator measuring itself."
    )
    chart(
        figures.known_null(cross),
        "known_null", "cross_fitting.csv")

    cols = st.columns(3)
    for col, arm in zip(cols, ("gc", "dn", "neg2")):
        frac = data.value(
            "cross_fitting.csv",
            f"2-mer fraction of the published value removed by cross-fitting, {arm} arm")
        pub = data.value("cross_fitting.csv", f"2-mer contribution as published, {arm} arm")
        col.metric(theme.PROTOCOL_LABEL[arm], f"{frac * 100:.1f}%",
                   help=f"Share of the two-stage value removed by cross-fitting. Two-stage "
                        f"reported {pub:.5f} where the truth is 0.")
    st.markdown(
        '<p class="src">the dinucleotide-matched arm exceeds 100% because cross-fitting carries '
        "the estimate slightly below zero, which is within its interval</p>",
        unsafe_allow_html=True)

    positives = [int(data.value("estimator_floor.csv",
                                f"datasets with a positive floor, {a} arm"))
                 for a in ("gc", "dn", "neg2")]
    box("note",
        "This is the argument for the cross-fitted estimator. On a question whose answer is "
        f"zero, the conventional estimator returns a positive value in {positives[0]}, "
        f"{positives[1]} and {positives[2]} of 94 datasets across the three protocols.")
    src("estimator_floor.csv")

    heading("The same comparison on the real signal", "pub")
    st.markdown(
        "The 4-mer does carry information the baseline lacks, so here cross-fitting does not "
        "collapse the estimate. It moves it slightly **up**."
    )
    chart(
        figures.animated_estimator_walk(cross),
        "animated_estimator_walk", "cross_fitting.csv")
    with st.expander("The same two estimators, side by side and still"):
        chart(
            figures.estimator_comparison(cross),
            "estimator_comparison", "cross_fitting.csv")
    st.markdown(
        "Which is why the cross-fitted value is the primary estimand rather than a correction to "
        "a too-large number. The span narrows from "
        f"{data.value('cross_fitting.csv', '4-mer three-arm span, as published'):.2f}x to "
        f"{data.value('cross_fitting.csv', '4-mer three-arm span, fully cross-fitted'):.2f}x, "
        "and the protocol dependence survives both."
    )

    with st.expander("Estimator floor table"):
        st.dataframe(data.table("estimator_floor.csv"), use_container_width=True,
                     hide_index=True)
        src("estimator_floor.csv")


def view_external(flt) -> None:
    page("External validation")
    st.markdown(
        "A held-out benchmark from a separate published deposit, sharing no dataset with the 94 "
        "analysed here."
    )
    if not guard("external_replication.csv", "external_replication_per_dataset.csv"):
        return

    external = data.table("external_replication.csv")
    span, lo, hi = data.interval("external_replication.csv",
                                 "SPAN across their two negative-set constructions")

    c1, c2, c3 = st.columns(3)
    c1.metric("Datasets", f"{int(data.value('external_replication.csv', 'datasets'))}",
              help=term("held-out benchmark"))
    c2.metric("Proteins", f"{int(data.value('external_replication.csv', 'proteins'))}",
              help="The resampled unit for the bootstrap.")
    c3.metric("Dataset overlap", "0", help="Zero by construction against the 94-dataset panel.")

    heading("What replicated", "pub")
    box("good",
        f"<strong>Protocol dependence replicated.</strong> Span {span:.2f}x (95% CI {lo:.2f} to "
        f"{hi:.2f}), against pre-fixed criteria of span above 1.5 and interval lower bound above "
        "1.2. Excluding every protein shared with this study's panel the span is "
        f"{data.value('external_replication.csv', _NO_SHARED_SPAN):.2f}x and still meets them.")
    chart(
        figures.external_arms(data.table("external_replication_per_dataset.csv")),
        "external_arms", "external_replication_per_dataset.csv", zoom=True)

    heading("What did not replicate", "pub")
    box("bad",
        "<strong>The inverse relation between apparent difficulty and contribution did not "
        "replicate.</strong> On this benchmark the arm with the higher composition baseline also "
        "has the higher contribution, so the two move in the <em>same</em> direction. Within the "
        "94-dataset panel they move in opposite directions.")
    left, right = st.columns(2)
    with left:
        st.plotly_chart(figures.external_levels(external, "baseline"),
                        use_container_width=True, theme=None, config=CHART)
    with right:
        st.plotly_chart(figures.external_levels(external, "contribution"),
                        use_container_width=True, theme=None, config=CHART)
    st.markdown(f'<p class="plain">{copy.PLAIN["external_levels"]}</p>',
                unsafe_allow_html=True)
    src("external_replication.csv")
    box("note",
        "Both charts read the same way: negative-1 is higher on the left and higher on the "
        "right. Two findings did not transport. The pooled baseline-to-contribution relation is "
        "null at their sample size, and the direction of the level relation is reversed. The "
        "manuscript does not present the reversal as a confirmed result, because the ordering "
        "was not predicted in advance.")

    with st.expander("Full external replication table"):
        st.dataframe(external, use_container_width=True, hide_index=True)
        src("external_replication.csv")


def view_reproducibility(flt) -> None:
    page("Reproducibility")
    st.markdown("How each committed table was produced, and what the release does not claim.")
    if not guard("PROVENANCE.csv"):
        return

    provenance = data.table("PROVENANCE.csv")
    box("warn", "<strong>Verification, not end-to-end reproduction.</strong> " + copy.REPRO_CAVEAT)

    heading("Provenance of the committed tables", "pub")
    chart(
        figures.provenance_mix(provenance),
        "provenance_mix", "PROVENANCE.csv")

    st.markdown("##### What the classes mean")
    st.markdown(
        "- **raw-reproducible** — regenerated from raw inputs by a `run.sh` stage.\n"
        "- **evidence-recomputable** — regenerated from the committed evidence, not raw reads.\n"
        "- **frozen-cache** — reproduced from a committed cached intermediate.\n"
        "- **frozen-only** — committed output; its script exists but is not wired into `run.sh`.\n"
        "- **cloud-produced** — produced by a cloud job; reproducing it needs that job rerun.\n"
        "- **unattributed** — no producing script recorded. Counted rather than hidden."
    )
    box("note",
        "Each table carries a sha256 in <code>PROVENANCE.csv</code>, so a changed table is "
        "detectable. That is an integrity check on the artefact, which is a different and weaker "
        "claim than reproducing it from raw data.")

    heading("Filter the provenance table")
    classes = sorted(provenance["status"].dropna().unique())
    picked = st.multiselect("Provenance class", classes, default=classes)
    shown = provenance[provenance["status"].isin(picked)]
    st.markdown(f'<p class="src">{len(shown)} of {len(provenance)} tables</p>',
                unsafe_allow_html=True)
    st.dataframe(shown[["table", "producing_script", "status", "bytes"]],
                 use_container_width=True, hide_index=True, height=340)
    src("PROVENANCE.csv")

    heading("Architecture")
    st.markdown(
        "The implemented execution and evidence paths, committed at "
        "`docs/assets/rbp-architecture-overview.png`. Rendered unmodified."
    )
    asset = data.REPO_ROOT / "docs" / "assets" / "rbp-architecture-overview.png"
    if asset.exists():
        st.image(str(asset), use_container_width=True)
        st.markdown(
            '<p class="src">docs/assets/rbp-architecture-overview.png  ·  eight further '
            "diagrams in docs/architecture.md</p>", unsafe_allow_html=True)
    else:
        st.info("Architecture asset not found.")


# -------------------------------------------------------------------------------------- shell

GUIDE_KEY = {"Overview": "Overview"}

VIEWS = {
    "Overview": view_overview,
    "Why this exists": view_why,
    "How negatives are built": view_construction,
    "Protocol sensitivity": view_protocol,
    "Model comparison": view_models,
    "Dataset explorer": view_explorer,
    "Cross-fitting": view_crossfitting,
    "Does the fix work?": view_recommendation,
    "External validation": view_external,
    "Reproducibility": view_reproducibility,
}


def main() -> None:
    with st.sidebar:
        st.markdown(
            '<div class="brand">'
            '<div class="brand-k">RBP</div>'
            '<div class="brand-t">Benchmark<br>Calibration</div>'
            "</div>",
            unsafe_allow_html=True,
        )
        st.markdown(theme.RULE, unsafe_allow_html=True)
        # One radio, so selection stays single and Streamlit keeps the state. The numbering and
        # the section headers are CSS on nth-of-type; interleaving real markdown between options
        # is not possible with a single widget, and splitting it into several radios would allow
        # two views to be selected at once.
        # ?view=N deep-links a view. Worth having on its own so a link can point at the
        # cross-fitting argument rather than at the front door, and it is the only way to
        # screenshot a specific view headlessly, which is how the white chart backgrounds and
        # two false claims survived as long as they did.
        names = list(VIEWS)
        requested = st.query_params.get("view")
        try:
            default = max(0, min(len(names) - 1, int(requested) - 1))
        except (TypeError, ValueError):
            default = 0
        choice = st.radio(
            "View", names, index=default, label_visibility="collapsed",
            format_func=lambda name: f"{names.index(name) + 1:02d}   {name}",
        )
        st.markdown("---")
        uses_filter = copy.GUIDE[GUIDE_KEY.get(choice, choice)]["filters"] != "static"
        try:
            flt = sidebar_filters(data.table("three_arm_per_dataset.csv"), "n_dn", uses_filter)
        except data.SchemaError:
            flt = NO_FILTER
        state = copy.FILTER_STATE[copy.GUIDE[GUIDE_KEY.get(choice, choice)]["filters"]][0]
        st.markdown(
            f'<div class="navfoot">{list(VIEWS).index(choice) + 1} of {len(VIEWS)}'
            f'<span>{state or ""}</span></div>',
            unsafe_allow_html=True,
        )
        st.markdown("---")
        with st.expander("Glossary"):
            for name, definition in copy.GLOSSARY.items():
                st.markdown(f"**{name}** — {definition}")
        st.markdown(
            '<p class="src">reads committed CSV only · never writes<br>'
            f'<a href="{PAPER_DOI}">preprint</a> · <a href="{CODE_DOI}">code</a></p>',
            unsafe_allow_html=True)
    VIEWS[choice](flt)


main()
