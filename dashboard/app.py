"""RBP Benchmark Calibration Explorer.

A reader over the committed tables in `results/tables/`. It opens CSV files, filters them in
memory and draws them. It does not run an analysis, recompute an estimate, reach any network or
write to the repository.

Run it with `streamlit run dashboard/app.py` from the repository root.
"""

from __future__ import annotations

import streamlit as st
from rbp_dashboard import copy, data, figures, theme

st.set_page_config(
    page_title="RBP Benchmark Calibration Explorer",
    page_icon=None,
    layout="wide",
    initial_sidebar_state="expanded",
)
st.markdown(theme.CSS, unsafe_allow_html=True)

_NO_SHARED_PROTEIN_SPAN = "SPAN excluding every protein shared with our panel"

PAPER_DOI = "https://doi.org/10.20944/preprints202609.0883.v1"
CODE_DOI = "https://doi.org/10.5281/zenodo.22679284"


def src(*names: str) -> None:
    """Name the file every figure came from. A chart without its source is an assertion."""
    listed = "  ·  ".join(f"results/tables/{n}" for n in names)
    st.markdown(f'<p class="src">source: {listed}</p>', unsafe_allow_html=True)


def note(text: str) -> None:
    st.markdown(f'<div class="note">{text}</div>', unsafe_allow_html=True)


def warn(text: str) -> None:
    st.markdown(f'<div class="warn">{text}</div>', unsafe_allow_html=True)


def term(name: str) -> str:
    return copy.GLOSSARY[name]


def guard(*tables: str) -> bool:
    """Validate before drawing. A missing column should say so, not raise inside a figure."""
    try:
        for name in tables:
            data.load(name)
    except data.SchemaError as exc:
        st.error(f"This view cannot render: {exc}")
        return False
    return True


# ----------------------------------------------------------------------------- views

def view_overview() -> None:
    st.title("RBP Benchmark Calibration Explorer")
    st.markdown(
        "How much does an RNA sequence model add beyond plain nucleotide composition, and how "
        "much does that answer depend on how the *unbound* windows were chosen?"
    )
    if not guard("cross_fitting.csv", "three_arm_contrast.csv", "panel_summary.csv"):
        return

    panel = data.table("panel_summary.csv")
    span, lo, hi = data.interval("cross_fitting.csv", "4-mer three-arm span, fully cross-fitted")
    two_stage = data.value("cross_fitting.csv", "4-mer three-arm span, as published")

    removed = min(
        data.value("cross_fitting.csv",
                   f"2-mer fraction of the published value removed by cross-fitting, {arm} arm")
        for arm in ("gc", "dn", "neg2")
    )

    c1, c2, c3, c4 = st.columns(4)
    c1.metric(
        "Primary span", f"{span:.2f}x",
        help="Ratio of the largest to the smallest nested contribution across the three "
             "negative-set protocols, for the 4-mer logistic regression, under the "
             "cross-fitted estimator. " + term("span"),
    )
    c2.metric(
        "Datasets", f"{int((panel['in_both_arms'] == True).sum())}",  # noqa: E712
        help=copy.PANEL_CAVEAT,
    )
    c3.metric("Protocols", "3", help="GC-matched, dinucleotide-matched and bias-aware.")
    c4.metric(
        "Known-null bias removed", f">={removed * 100:.0f}%",
        help="On a constructed test whose true answer is zero, cross-fitting removes at least "
             "this share of what the conventional estimator reports. Smallest of the three "
             "protocols.",
    )
    st.markdown(
        f"<p class='src'>95% CI {lo:.2f} to {hi:.2f}, protein-clustered bootstrap  ·  "
        f"two-stage estimator gives {two_stage:.2f}x on the same data</p>",
        unsafe_allow_html=True,
    )

    st.markdown("## The finding")
    st.markdown(
        "Hold the model class, the source peaks, the chromosome-blocked folds and the estimator "
        "fixed. Change how the negative windows are built. The measured contribution moves "
        f"**{span:.2f}-fold**, and it moves *against* the apparent score: the protocol that "
        "looks easiest by headline AUROC is the one where the model contributes most."
    )
    st.plotly_chart(
        figures.protocol_contribution(data.table("three_arm_contrast.csv")),
        use_container_width=True, config={"displayModeBar": False},
    )
    src("three_arm_contrast.csv")

    note(copy.NOT_ONLY_NEGATIVES)
    warn(
        "<strong>What this dashboard is.</strong> " + copy.NOT_MONITORING
        + " Every number shown is read from a committed table; none is recomputed here."
    )

    st.markdown("## Scope")
    left, right = st.columns(2)
    with left:
        st.markdown(
            f"- **Paper:** [Preprints.org]({PAPER_DOI}) (preprint; not yet peer reviewed)\n"
            f"- **Code and evidence:** [Zenodo]({CODE_DOI})\n"
            "- **Panel:** 94 ENCODE eCLIP datasets carrying all three protocols\n"
            "- **Resampling:** protein-clustered bootstrap"
        )
    with right:
        st.markdown(
            "- **Primary estimator:** cross-fitted\n"
            "- **Comparability estimator:** two-stage, what the surveyed literature computes\n"
            "- **Held-out benchmark:** 135 datasets, no dataset shared\n"
            "- **Model classes:** 4-mer, CNN, SpliceBERT"
        )


def view_protocol() -> None:
    st.title("Protocol sensitivity")
    st.markdown(
        "Three ways of choosing the unbound windows, applied to the same 94 datasets with "
        "everything else held fixed."
    )
    if not guard("three_arm_contrast.csv", "three_arm_per_dataset.csv"):
        return

    contrast = data.table("three_arm_contrast.csv")
    per_dataset = data.table("three_arm_per_dataset.csv")

    with st.expander("What the three protocols are"):
        for key in ("GC-matched", "dinucleotide-matched", "bias-aware"):
            st.markdown(f"**{key}.** {term(key)}")

    st.markdown("## Apparent AUROC and its baseline")
    st.markdown(
        "Two quantities on one scale. The dark bar is what plain composition scores on its own; "
        "the blue bar is composition plus the model's score, which is what a paper would print."
    )
    st.plotly_chart(figures.protocol_levels(contrast), use_container_width=True,
                    config={"displayModeBar": False})
    src("three_arm_contrast.csv")

    st.markdown("## What the model actually adds")
    st.markdown(
        "The difference between those two bars, on its own axis. "
        f"*{term('nested contribution')}*"
    )
    st.plotly_chart(figures.protocol_contribution(contrast), use_container_width=True,
                    config={"displayModeBar": False})
    src("three_arm_contrast.csv")
    note(
        "Read the two charts together. The bias-aware protocol has the <em>highest</em> apparent "
        "AUROC and the <em>lowest</em> contribution. The dinucleotide-matched protocol is the "
        "reverse. A headline AUROC alone cannot distinguish them."
    )

    st.markdown("## Per dataset, not just the panel mean")
    st.markdown("Every dataset as a point. The panel means above are the summary of these.")
    st.plotly_chart(figures.contribution_distribution(per_dataset), use_container_width=True,
                    config={"displayModeBar": False})
    src("three_arm_per_dataset.csv")

    st.markdown("## Paired within dataset")
    left, right = st.columns([1, 2])
    with left:
        pair = st.selectbox(
            "Compare", ["dn vs gc", "neg2 vs gc", "dn vs neg2"], index=0,
            help="The same datasets appear in both arms, so the change is paired.",
        )
    mapping = {"dn vs gc": ("gc", "dn"), "neg2 vs gc": ("gc", "neg2"), "dn vs neg2": ("neg2", "dn")}
    a, b = mapping[pair]
    st.plotly_chart(figures.paired_slopes(per_dataset, a, b), use_container_width=True,
                    config={"displayModeBar": False})
    src("three_arm_per_dataset.csv")

    st.markdown("## The inverse relation")
    st.markdown(
        "Apparent AUROC on the horizontal axis, contribution on the vertical. Within this panel "
        "the two move in opposite directions. **This is the finding that did not replicate on the "
        "held-out benchmark** (see External validation)."
    )
    st.plotly_chart(figures.apparent_versus_contribution(per_dataset), use_container_width=True,
                    config={"displayModeBar": False})
    src("three_arm_per_dataset.csv")

    with st.expander("Numbers behind these charts"):
        st.dataframe(
            contrast[["check", "value", "ci_low", "ci_high", "n"]],
            use_container_width=True, hide_index=True,
        )
        src("three_arm_contrast.csv")


def view_models() -> None:
    st.title("Model comparison")
    st.markdown(
        "Three model classes, each measured under all three protocols on the same 94 datasets."
    )
    if not guard("three_arm_models.csv", "three_arm_models_per_dataset.csv"):
        return

    models = data.table("three_arm_models.csv")

    warn("<strong>Estimator notice.</strong> " + copy.ESTIMATOR_CAVEAT)

    st.markdown("## How much each model's answer moves with the protocol")
    st.plotly_chart(figures.model_spans(models), use_container_width=True,
                    config={"displayModeBar": False})
    src("three_arm_models.csv")
    st.markdown(
        "All three spans are far from 1.0, so protocol dependence is not an artefact of one "
        "model class. The spans are **not** comparable to each other as estimates of the same "
        "quantity: the 4-mer figure here is two-stage, shown so the three sit on one footing, "
        "while the paper's primary 4-mer span is the cross-fitted "
        f"{data.value('cross_fitting.csv', '4-mer three-arm span, fully cross-fitted'):.2f}x."
    )

    st.markdown("## Contribution within each protocol")
    st.plotly_chart(figures.model_by_protocol(models), use_container_width=True,
                    config={"displayModeBar": False})
    src("three_arm_models.csv")
    note(
        "All three model classes put their lowest contribution in the bias-aware arm. "
        "SpliceBERT contributes more than the others in every protocol, and its span is the "
        "narrowest, but it is still a "
        f"{data.value('three_arm_models.csv', 'splicebert three-protocol span'):.2f}-fold range."
    )

    with st.expander("Per-dataset values"):
        st.dataframe(data.table("three_arm_models_per_dataset.csv"),
                     use_container_width=True, hide_index=True)
        src("three_arm_models_per_dataset.csv")


def view_crossfitting() -> None:
    st.title("Cross-fitting")
    st.markdown("What happens when the estimator is tested on a question whose answer is known.")
    if not guard("cross_fitting.csv", "estimator_floor.csv"):
        return

    cross = data.table("cross_fitting.csv")

    with st.expander("The two estimators, in plain language", expanded=True):
        st.markdown(f"**Two-stage.** {term('two-stage estimator')}")
        st.markdown(f"**Cross-fitted.** {term('cross-fitted estimator')}")
        st.markdown(f"**Estimator floor.** {term('estimator floor')}")

    st.markdown("## The known-answer test")
    st.markdown(
        "Score each window with a **2-mer** model. The composition baseline already contains "
        "every 2-mer frequency, so a 2-mer model can add nothing: the true contribution is "
        "**zero by construction**. Anything an estimator reports here is measuring itself."
    )
    st.plotly_chart(figures.known_null(cross), use_container_width=True,
                    config={"displayModeBar": False})
    src("cross_fitting.csv")

    cols = st.columns(3)
    for col, arm in zip(cols, ("gc", "dn", "neg2")):
        frac = data.value(
            "cross_fitting.csv",
            f"2-mer fraction of the published value removed by cross-fitting, {arm} arm")
        pub = data.value("cross_fitting.csv", f"2-mer contribution as published, {arm} arm")
        col.metric(
            theme.PROTOCOL_LABEL[arm], f"{frac * 100:.1f}%",
            help=f"Share of the two-stage value removed by cross-fitting. The two-stage "
                 f"estimator reported {pub:.5f} where the truth is 0.",
        )
    st.markdown(
        "<p class='src'>the dinucleotide-matched arm exceeds 100% because cross-fitting carries "
        "the estimate slightly below zero, which is within its interval</p>",
        unsafe_allow_html=True,
    )
    positives = [
        int(data.value("estimator_floor.csv", f"datasets with a positive floor, {a} arm"))
        for a in ("gc", "dn", "neg2")
    ]
    note(
        "This is the argument for the cross-fitted estimator. On a question whose answer is zero, "
        "the conventional estimator returns a positive value in "
        f"{positives[0]}, {positives[1]} and {positives[2]} of 94 datasets across the three "
        "protocols."
    )
    src("estimator_floor.csv")

    st.markdown("## The same comparison on the real signal")
    st.markdown(
        "The 4-mer model does carry information the composition baseline lacks, so here "
        "cross-fitting does not collapse the estimate. It moves it slightly **up**."
    )
    st.plotly_chart(figures.estimator_comparison(cross), use_container_width=True,
                    config={"displayModeBar": False})
    src("cross_fitting.csv")
    st.markdown(
        "Which is why the paper reports the cross-fitted value as its primary estimand rather "
        "than as a correction to a too-large number. The span narrows from "
        f"{data.value('cross_fitting.csv', '4-mer three-arm span, as published'):.2f}x to "
        f"{data.value('cross_fitting.csv', '4-mer three-arm span, fully cross-fitted'):.2f}x, "
        "and the protocol dependence survives both."
    )

    with st.expander("Estimator floor table"):
        st.dataframe(data.table("estimator_floor.csv"),
                     use_container_width=True, hide_index=True)
        src("estimator_floor.csv")


def view_external() -> None:
    st.title("External validation")
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
    c3.metric("Dataset overlap", "0",
              help="Zero by construction against the 94-dataset three-arm panel.")

    st.markdown("## What replicated")
    st.success(
        f"**Protocol dependence replicated.** Span {span:.2f}x (95% CI {lo:.2f} to {hi:.2f}) "
        "across their two negative-set constructions, against pre-fixed criteria of span above "
        "1.5 and interval lower bound above 1.2. Excluding every protein shared with this "
        "study's panel, the span is "
        f"{data.value('external_replication.csv', _NO_SHARED_PROTEIN_SPAN):.2f}x "
        "and still meets those criteria."
    )
    st.plotly_chart(figures.external_arms(data.table("external_replication_per_dataset.csv")),
                    use_container_width=True, config={"displayModeBar": False})
    src("external_replication_per_dataset.csv")

    st.markdown("## What did not replicate")
    st.error(
        "**The inverse relation between apparent difficulty and contribution did not replicate.** "
        "On this benchmark the arm with the higher composition baseline also has the higher "
        "contribution, so the two move in the *same* direction. Within the 94-dataset panel they "
        "move in opposite directions."
    )
    left, right = st.columns(2)
    with left:
        st.plotly_chart(figures.external_levels(external, "baseline"),
                        use_container_width=True, config={"displayModeBar": False})
    with right:
        st.plotly_chart(figures.external_levels(external, "contribution"),
                        use_container_width=True, config={"displayModeBar": False})
    src("external_replication.csv")
    note(
        "Both charts read the same way: negative-1 is higher on the left and higher on the right. "
        "Two findings did not transport. The pooled baseline-to-contribution relation is null at "
        "their sample size, and the direction of the level relation is reversed. The manuscript "
        "does not present the reversal as a confirmed result, because the ordering was not "
        "predicted in advance."
    )

    with st.expander("Full external replication table"):
        st.dataframe(external, use_container_width=True, hide_index=True)
        src("external_replication.csv")


def view_reproducibility() -> None:
    st.title("Reproducibility")
    st.markdown("How each committed table was produced, and what the release does and does "
                "not claim.")
    if not guard("PROVENANCE.csv"):
        return

    provenance = data.table("PROVENANCE.csv")

    warn("<strong>Verification, not end-to-end reproduction.</strong> " + copy.REPRO_CAVEAT)

    st.markdown("## Provenance of the committed tables")
    st.plotly_chart(figures.provenance_mix(provenance), use_container_width=True,
                    config={"displayModeBar": False})
    src("PROVENANCE.csv")

    st.markdown("### What the classes mean")
    st.markdown(
        "- **raw-reproducible** — regenerated from raw inputs by a `run.sh` stage.\n"
        "- **evidence-recomputable** — regenerated from the committed evidence, not raw reads.\n"
        "- **frozen-cache** — reproduced from a cached intermediate that is committed.\n"
        "- **frozen-only** — committed output; its generating script is present but not wired "
        "into `run.sh`.\n"
        "- **cloud-produced** — produced by a cloud job. Reproducing it needs that job rerun.\n"
        "- **unattributed** — no producing script recorded. Counted rather than hidden."
    )
    note(
        "Each table carries a sha256 in <code>PROVENANCE.csv</code>, so a changed table is "
        "detectable. That is an integrity check on the committed artefact, which is a different "
        "and weaker claim than reproducing the artefact from raw data."
    )

    st.markdown("## Architecture")
    st.markdown(
        "The implemented execution and evidence paths, as committed in "
        "`docs/assets/rbp-architecture-overview.png`. Rendered here unmodified."
    )
    asset = data.REPO_ROOT / "docs" / "assets" / "rbp-architecture-overview.png"
    if asset.exists():
        st.image(str(asset), use_container_width=True)
        st.markdown(
            "<p class='src'>source: docs/assets/rbp-architecture-overview.png  ·  "
            "eight further diagrams in docs/architecture.md</p>", unsafe_allow_html=True)
    else:
        st.info("Architecture asset not found at docs/assets/rbp-architecture-overview.png")

    with st.expander("Provenance table"):
        st.dataframe(provenance, use_container_width=True, hide_index=True)
        src("PROVENANCE.csv")


# ----------------------------------------------------------------------------- shell

VIEWS = {
    "Overview": view_overview,
    "Protocol sensitivity": view_protocol,
    "Model comparison": view_models,
    "Cross-fitting": view_crossfitting,
    "External validation": view_external,
    "Reproducibility": view_reproducibility,
}


def main() -> None:
    with st.sidebar:
        st.markdown("### RBP Calibration")
        choice = st.radio("View", list(VIEWS), label_visibility="collapsed")
        st.markdown("---")
        with st.expander("Glossary"):
            for name, definition in copy.GLOSSARY.items():
                st.markdown(f"**{name}** — {definition}")
        st.markdown("---")
        st.markdown(
            f"<p class='src'>Preprint: <a href='{PAPER_DOI}'>10.20944/preprints202609.0883.v1</a>"
            f"<br>Code: <a href='{CODE_DOI}'>10.5281/zenodo.22679284</a>"
            "<br>Reads committed CSV only. Never writes.</p>",
            unsafe_allow_html=True,
        )
    VIEWS[choice]()


main()
