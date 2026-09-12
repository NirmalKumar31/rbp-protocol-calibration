# RBP Benchmark Calibration Explorer

A read-only dashboard over the committed tables in `results/tables/`. It exists to make the
study's result legible without reading a 55-page paper, and to be a portfolio artefact.

## What it is not

It does not run analyses, fit models, recompute estimates, or touch any cloud service, and
it makes no network request for data. See **Fonts and network** below for the one request the
browser does make. It opens CSV files, filters them in memory and draws them. **It is not production
monitoring.** Every number it displays is read from a committed table.

## Run it

From the repository root, in a virtual environment that is **not** the study environment:

```
python -m venv .venv-dashboard
source .venv-dashboard/bin/activate
pip install -r dashboard/requirements.txt
streamlit run dashboard/app.py --server.address 127.0.0.1
```

`--server.address 127.0.0.1` binds it to this machine only. Streamlit otherwise listens on all interfaces, which puts it on your local network.

It opens on `http://localhost:8501`. The dependencies are pinned exactly and kept out of
`pyproject.toml` so that a Streamlit or Plotly upgrade can never alter a published number.

## Deploying it

The app is a read-only viewer over files already in this repository, so a host only needs the
repository and `dashboard/requirements.txt`. Nothing is written, no secret is required, and
there is no database.

**Streamlit Community Cloud**

| setting | value |
|---|---|
| repository | `NirmalKumar31/rbp-protocol-calibration` |
| branch | the branch carrying `dashboard/` |
| main file path | `dashboard/app.py` |
| Python version | 3.13 |

`requirements.txt` lives at `dashboard/requirements.txt`, **not** at the repository root. Point
the host at it explicitly if it does not find it; a host that silently installs nothing will
fail on `import streamlit` rather than on anything informative.

`.streamlit/config.toml` is at the **repository root**, not inside `dashboard/`. Streamlit reads
`$CWD/.streamlit/config.toml` and the app is launched from the root, so a config under
`dashboard/` is never read. It lived there unread for several releases: the dark surface comes
from injected CSS and looked correct, while Streamlit's own widgets kept their defaults. Run
`streamlit config show` to see what is actually in force.

### What deploying costs

**The app reads 18 CSV files totalling about 0.1 MB, plus one 550 KB image. Cloning this
repository fetches roughly 94 MB across 4,000 files.** Everything else is the scientific
evidence, the manuscript and its figures. That is not a failure, and Streamlit Community Cloud
allows it, but it means slow cold starts and it puts the full evidence tree on a third-party
host. The repository is already public and the preprint is posted, so nothing becomes visible
that was not already, but it is a decision rather than a detail.

If that matters, the alternative is a second repository holding only `dashboard/`, the 18
tables and the one image, at a few megabytes. The cost is that the tables would then be copies,
and a copy can drift from the release it was taken from, which is the defect class this whole
project exists to catch. **Deploying from this repository keeps the dashboard reading the same
files the paper cites.** That is the reason it is set up this way.

## Tests

```
bash dashboard/ci.sh
```

That runs ruff, the tests, and a render of all ten views at four filter widths including an
empty one. It is wired into `.github/workflows/ci.yml` as the **`dashboard`** job, which installs
`dashboard/requirements.txt` and runs this script. A separate job because the study environment
must not acquire Streamlit or Plotly: a Plotly bump cannot be allowed to change a published
number. The `test` and `full-suite` jobs run `pytest tests`, which never reaches
`dashboard/tests`, so before this job the dashboard suite ran nowhere.

Three things are gated:

1. **Schemas.** Every table a figure reads exists, is non-empty, and has the columns that figure
   needs. A missing column fails here rather than inside a chart.
2. **Headline values.** The figures displayed on the Overview are compared against the
   `check` rows in the committed tables. If a table changes and the dashboard's story does not,
   this fails.
3. **Read-only.** No module in `rbp_dashboard/` may contain a write call. Checked by parsing the
   source, not by convention.

## Views

| View | Source tables |
|---|---|
| Overview | `cross_fitting.csv`, `three_arm_contrast.csv`, `panel_summary.csv` |
| Dataset explorer | `three_arm_per_dataset.csv`, `three_arm_models_per_dataset.csv` |
| Protocol sensitivity | `three_arm_contrast.csv`, `three_arm_per_dataset.csv` |
| Model comparison | `three_arm_models.csv`, `three_arm_models_per_dataset.csv` |
| Cross-fitting | `cross_fitting.csv`, `estimator_floor.csv` |
| External validation | `external_replication.csv`, `external_replication_per_dataset.csv` |
| Reproducibility | `PROVENANCE.csv`, `docs/assets/rbp-architecture-overview.png` |

Every chart names its own source file beneath it, in the interface.

## Scientific labelling

The dashboard is constrained to state the study's claims at their actual strength:

- The **4.84x span is the primary, cross-fitted** result, and applies to the **4-mer logistic
  regression**. The two-stage figure for the same quantity is 5.42x and is shown as the
  comparability analysis, because it is what the surveyed literature computes.
- The **CNN and SpliceBERT spans are two-stage and exploratory.** Neither was cross-fitted. Every
  view that compares model classes carries that notice.
- The **95% bias reduction** is measured on the constructed **2-mer known-null test**, where the
  true contribution is zero by construction, across the three protocols. It is not a general
  claim about the estimator.
- **95 datasets** were selected for the panel; **94** carry all three protocols and are the
  denominator for three-protocol comparisons.
- Negative-set construction is what varies, but it is **not the only thing that changes**: models
  and baselines are refitted per protocol and retained positive subsets differ slightly.
- Reproducibility is described as **verification against committed evidence**, not complete
  raw-to-result or bit-for-bit reproduction.
- The held-out benchmark **replicated protocol dependence** and **did not replicate** the inverse
  relation between apparent difficulty and contribution. Both are stated.

## Interactivity

Sidebar filters on cell line, protein and dataset size narrow **every** per-dataset chart at
once. The Dataset explorer adds a sortable heatmap of all 94 datasets against the three
protocols, an adjustable ranking, and a drill-down into any single dataset across all three
model classes. The Model comparison view carries a genuine Plotly frame animation stepping
through the model classes, so the protocol ordering can be watched holding while the bars move.

### The rule that constrains it

Filtering is the interactive part. The published estimand does not move when you filter.

- A value read from a committed table is marked **PUBLISHED** and carries its
  protein-clustered bootstrap interval.
- A value computed live from a filtered subset is marked **LIVE / DESCRIPTIVE**, carries **no
  interval**, and is never called an estimate.

`data.descriptive()` returns only `n`, `mean`, `median`, `min` and `max`. A test asserts no
interval key is present, because if one ever were, some figure would eventually draw an error
bar on a filtered mean and present it as an estimate. Where the two appear together, the
published bar is solid and the filtered bar is hatched.

Single-dataset values are shown because they are in the committed table, with a notice that the
study's unit of inference is the panel and one dataset carries no interval at all.

## Colour

Dark surface, `#10141C`. Three data colours, measured rather than chosen: pairwise OKLab
separation under normal vision and under simulated protanopia, deuteranopia and tritanopia
(Machado 2009 matrices at full severity), plus contrast against the surface.

| | |
|---|---|
| worst normal-vision separation | dE 18.8 (floor 15) |
| worst CVD separation | dE 17.6 (target 8) |
| contrast against surface | 3.26:1 to 10.54:1 |

The colours are steel `#A8C8E0`, brass `#BE9752` and slate `#55697F`.

**Separation is carried by lightness, not saturation, and that is why the palette looks calm.**
Four desaturated sets were measured first and all four failed: muting a hue removes chroma, and
chroma is what CVD separation depends on, so a soft blue and a soft clay collapse into each
other under deuteranopia. Spreading three colours across a wide lightness range buys separation
that survives every simulation while keeping every colour muted.

Two earlier palettes are worth recording because both were wrong in instructive ways. A light
theme used blue, amber and slate on white; **a palette validated on white does not transfer to
near-black**, so the dark theme was re-measured from scratch rather than inverted. A first dark
attempt used saturated cyan, amber and rose, which measured well and read as neon.

`theme.NUCLEOTIDE` (A, C, G, U) is the header motif and the background field, and is
**decorative only**. Those four fail CVD separation against each other, which is acceptable
where nothing depends on telling them apart and would be a defect in a legend. A test asserts
`figures.py` never references it.

## Fonts and network

The dashboard reads only local CSV files and makes no network request for data. It does load
**Inter and JetBrains Mono from Google Fonts** over the network, so the browser contacts
`fonts.googleapis.com` on first paint. Nothing scientific leaves the machine; if that request is
unacceptable, delete the `@import` in `rbp_dashboard/theme.py` and the stack falls back to the
system sans and monospace faces already listed in `FONT` and `MONO`.
