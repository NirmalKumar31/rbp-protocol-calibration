# What a sequence model appears to contribute is set by how the negatives were built

A calibration study across 94 paired ENCODE eCLIP datasets. We hold the model, the peak set, the
chromosome-blocked folds and the estimator fixed, change only how negative windows are
constructed, and measure each model's **nested contribution**: the out-of-fold AUROC of a
logistic model on 19 composition features plus the model's score, minus the AUROC of those
features alone.

The number moves **5.42-fold** (95% CI 4.43 to 6.58) for a 4-mer logistic regression, and 3.7 to
7.4-fold across three model classes on identical rows. A model's apparent AUROC moves the
opposite way. And the estimator itself returns **+0.011 to +0.014** when the true contribution is
zero by construction, which is 91% of the smallest arm's reported value.

> **Report the composition-only AUROC obtained under the same protocol alongside every headline
> AUROC. Do not compare contributions measured under different protocols.**

## Check it in thirty seconds, offline

No cloud account, no credentials, no data download. Everything the paper asserts is checked
against committed tables.

```bash
git clone <this repo> && cd rbp-repro
python -m pip install -e .
PYTHONPATH=src python scripts/verify.py --local results/tables   # 937/937
PYTHONPATH=src python -m pytest tests -q                          # 699 passed
```

`verify.py` re-derives every published value from the committed result tables and fails if any
disagrees with `config/golden.yaml`. Two of its assertions are end-to-end rebuilds rather than
comparisons against a record: 285 AUROCs recomputed from committed per-window scores, and the
headline contrast recomputed from raw sequence.

## What the study found

| | |
|---|---|
| **The protocol moves the measurement** | Nested contribution for one 4-mer: **+0.0663** dinucleotide-matched, **+0.0265** GC-matched, **+0.0122** bias-aware (negatives are other RBPs' sites). Apparent AUROC moves the other way, 0.798 to 0.689, in 94 of 94 datasets |
| **It holds for three model classes** | Spans of 5.42, 7.42 and 3.72 for a 4-mer, a 7089-parameter CNN and a fine-tuned SpliceBERT, on identical rows and folds |
| **It is not an AUROC artefact** | The ordering holds on five estimands including unbounded deviance. The magnitude is scale-specific: 5.42-fold in AUROC, about 2.1-fold on unbounded scales |
| **Shuffling removes the baseline entirely** | Dinucleotide-shuffled negatives pin the composition baseline at exactly **0.5000** on all 94 datasets, so the contribution becomes the model's own AUROC less a half. Across four constructions the span is 20.62-fold |
| **The estimator has a floor** | Applied to a model whose information the baseline already contains, so the truth is zero, it returns **+0.0119 / +0.0137 / +0.0111**. Nearly flat across arms, so the span survives; but 90.6% of the bias-aware arm's value, so that level does not |
| **The baseline's order matters too** | Raising it to order three removes most of a 4-mer's contribution and a third of SpliceBERT's; at order four the baseline overfits and the estimator's error exceeds most published increments |
| **Nobody reports the baseline** | Of seven widely used methods and benchmarks, five build negatives by relocating genomic intervals, which leaves composition unconstrained, and **none** reports a composition-only AUROC |

## Reproduce it in full

The offline check above needs nothing. Regenerating the result tables from raw data needs the
window store and, for the neural arms, a GPU.

```bash
export GOOGLE_CLOUD_PROJECT=your-project
./run.sh preflight        # spends nothing, gates everything
./run.sh all              # pauses before every paid stage
```

Full procedure: **[docs/REPRODUCE.md](docs/REPRODUCE.md)**.
Why dataset counts differ between analyses: **[docs/PANELS.md](docs/PANELS.md)**.

Cost: about **$5** of GCP and **$19** of Modal GPU for the neural sweeps, which are 95% of it.
Every paid stage asks first.

## Layout

```
manuscript/     the paper and its figures
scripts/        one analysis per file; each writes a table under results/tables/
src/rbp/        the library the scripts import
tests/          699 tests, no network or cloud required
config/         params.yaml (the study's settings), golden.yaml (expected values)
results/tables/ every number in the paper
data/evidence/  per-window out-of-fold scores for all three model classes
```

## Design rules

1. **Verification is a stage**, not an afterthought. Reproducibility that is not checked is not
   reproducibility.
2. **No hardcoded project id.** Everything resolves through `rbp.utils.cloud`, and a test fails
   the build if a literal reappears.
3. **The panel is an artefact, not a flag** (`manifest/study_panel.tsv`), written once.
4. **Task counts come from manifests**, never typed by hand.
5. **Completion markers are written last**, so an interrupted stage redoes its work rather than
   being skipped.

## Licence and citation

Code and derived data under MIT and CC BY 4.0 respectively; see `LICENSE`. Intermediate window
tables containing genomic sequence are not redistributed and are regenerated from the ENCODE
accessions in Supplementary Table S1.
