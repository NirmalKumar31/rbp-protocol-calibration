# RBP Benchmark Calibration Explorer

A read-only dashboard over the committed tables in `results/tables/`. It exists to make the
study's result legible without reading a 55-page paper, and to be a portfolio artefact.

## What it is not

It does not run analyses, fit models, recompute estimates, reach a network, or touch any cloud
service. It opens CSV files, filters them in memory and draws them. **It is not production
monitoring.** Every number it displays is read from a committed table.

## Run it

From the repository root, in a virtual environment that is **not** the study environment:

```
python -m venv .venv-dashboard
source .venv-dashboard/bin/activate
pip install -r dashboard/requirements.txt
streamlit run dashboard/app.py
```

It opens on `http://localhost:8501`. The dependencies are pinned exactly and kept out of
`pyproject.toml` so that a Streamlit or Plotly upgrade can never alter a published number.

## Tests

```
python -m pytest dashboard/tests -q
python -m ruff check dashboard
```

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

## Colour

Three categorical colours, checked rather than chosen by eye: pairwise OKLab separation under
normal vision and under simulated protanopia, deuteranopia and tritanopia (Machado 2009 matrices
at full severity), plus contrast against the surface.

| | |
|---|---|
| worst normal-vision separation | dE 20.4 (floor 15) |
| worst CVD separation | dE 14.4 (target 8) |
| weakest contrast against surface | 3.19:1 |

An earlier blue/amber/teal set failed at dE 6.4 because blue and teal collapse together under
tritanopia. The surface is pinned light in `.streamlit/config.toml` because that is the surface
the palette was validated against.
