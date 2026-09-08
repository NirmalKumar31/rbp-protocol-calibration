# Apparent sequence-model contribution depends strongly on negative-set construction

A calibration study across 94 paired ENCODE eCLIP datasets. We hold the model class and its
hyperparameters, the peak set, the chromosome-to-fold map and the estimator implementation
fixed, change how negative windows are constructed, and measure each model's **nested
contribution**: the out-of-fold AUROC of a logistic model on 19 composition features plus the
model's score, minus the AUROC of those features alone.

Two things are not held fixed and the abstract says so: model and baseline are **refitted per
protocol**, and because each matcher rejects the positives it cannot pair, the retained
positives differ very slightly between arms. Restricting to their intersection drops 0.23% of
positives and moves the contrast +0.0398 to +0.0401.

The number moves **4.84-fold** (95% CI 3.98 to 5.81) for a 4-mer logistic regression under the
cross-fitted estimator this paper recommends, and **5.42-fold** (4.43 to 6.58) under the
two-stage estimator the literature actually computes. Both are reported throughout; the
cross-fitted one is primary. Across three model classes the two-stage span runs 3.7 to 7.4-fold,
on identical rows within each arm. A model's apparent AUROC moves the opposite way.

Why there are two. The two-stage estimator returns **+0.011 to +0.014** when the true
contribution is zero by construction, which is 90.4% of the smallest arm's reported value.
Cross-fitting the score covariate cuts that floor by at least 95% and returns the zero the
construction requires. The effect survives the correction, and the neural spans are two-stage
only, so they stay exploratory.

> **Report the composition-only AUROC obtained under the same protocol alongside every headline
> AUROC. Do not compare contributions measured under different protocols.**

## If you are reviewing this

The primary claim is that a model's **measured** contribution over a composition baseline
depends strongly on how the negative windows were built: 4.84-fold cross-fitted across three
protocols for a 4-mer, 5.42-fold under the estimator the literature uses, while its apparent
AUROC moves the opposite way. A two-way decomposition of the nine
train-by-evaluate combinations attributes most of the movement to the **evaluation** protocol
rather than the fitted model. Two decompositions, which are different estimands and not two
weightings of one: averaging each dataset's own normalised shares gives 63% (CI 57 to 68)
against 15% for training, and decomposing the matrix of panel means gives 81% (73 to 87)
against 9%. Leave-one-protein-out moves the larger share by at most 3.0 points, so no one
protein carries it. Most of what the protocol moves is the measurement (`sec:transport`).
It is descriptive, not causal.

- **Verifiable offline, in one command, in under a minute:** every published number, against
  committed tables. That is regression checking, not independent reproduction.
- **Not regenerable without cloud and raw data:** the neural sweeps for two of three protocols,
  and anything needing the 2.9 GB window store.
  `results/tables/PROVENANCE.csv` says which of the two every released table is.
- **Known to be incomplete:** the cross-fitted estimator we recommend is computed for the k-mer
  classes only; the two neural spans come from the estimator we say to replace, and are
  exploratory. Intervals condition on one negative draw and one fold partition.

The panel is 95 datasets; 94 carry all three protocols, and the one that does not
(NCBP2 in K562) is named in Supplementary Table S1. Both counts are correct and appear
throughout for different quantities.

## Check it in thirty seconds, offline

No cloud account, no credentials, no data download. 1153 numeric assertions are checked
against committed tables, of which **1015 belong to this paper** and 136 to an earlier
variant-scoring study whose code and evidence are still here and still pass. The verifier prints
that split on every run, because one total covering two papers is not this paper's evidence. That is a regression gate on the published values, not a proof that
each is attached to the right claim; the Limitations section says what it does not cover.

```bash
git clone https://github.com/NirmalKumar31/rbp-protocol-calibration.git && cd rbp-protocol-calibration
python -m pip install -e . -c constraints.txt   # no torch: the neural stack is an extra

PYTHONPATH=src python scripts/verify.py --local results/tables   # 1153/1153
PYTHONPATH=src python -m pytest tests -q \
  --ignore=tests/unit/test_models.py --ignore=tests/unit/test_train_folds.py

# Everything CI runs, in the order that works. WITHOUT torch this prints PREFLIGHT PARTIAL
# and names the steps it could not run; it will not report CLEAN over a subset.
./scripts/preflight.sh          # check only
./scripts/preflight.sh --fix    # also sync derived counts and refresh the manifests

# For the FULL release gate, and to reproduce the sweeps, install the neural extra first:
python -m pip install -e '.[neural]' -c constraints.txt
./scripts/preflight.sh          # now PREFLIGHT CLEAN means the whole gate
```

`verify.py` re-derives every published value from the committed result tables and fails if any
disagrees with `config/golden.yaml`. Two of its assertions are end-to-end rebuilds rather than
comparisons against a record: 285 AUROCs recomputed from committed per-window scores, and the
headline contrast recomputed from raw sequence.

## What the study found

| | |
|---|---|
| **The protocol moves the measurement** | Nested contribution for one 4-mer: **+0.0663** dinucleotide-matched, **+0.0265** GC-matched, **+0.0122** bias-aware (negatives are other RBPs' sites). Apparent AUROC moves the other way, 0.798 to 0.688, in 94 of 94 datasets |
| **It holds for three model classes** | Spans of 5.42, 7.42 and 3.72 for a 4-mer, a 7089-parameter CNN and a fine-tuned SpliceBERT. Within each arm all three are scored on identical rows and folds; rows differ between arms by construction |
| **It is not an AUROC artefact** | The ordering holds on five estimands including unbounded deviance. The magnitude is scale-specific: 5.42-fold in AUROC, about 2.1-fold on unbounded scales |
| **Shuffling removes the baseline entirely** | Dinucleotide-shuffled negatives pin the composition baseline at exactly **0.5000** on all 94 datasets, so the contribution becomes the model's own AUROC less a half. Across four constructions the span is 20.62-fold |
| **The estimator has a floor** | Applied to a model whose information the baseline already contains, so the truth is zero, it returns **+0.0119 / +0.0137 / +0.0111**. Nearly flat across arms, so the span survives; but 90.4% of the bias-aware arm's value, so that level does not |
| **The floor is removable** | It is the outer-fold route, not conditioning. Cross-fitting the covariate cuts it by **at least 95%** and lands within 5e-4 of the known zero. The span goes 5.42 to **4.84**. Measured for the k-mer classes; the CNN and SpliceBERT would need about $115, twice one sweep of both models across the three arms |
| **The baseline's order matters too** | Raising it to order three removes most of a 4-mer's contribution and a third of SpliceBERT's; at order four the baseline overfits and the estimator's error exceeds most published increments |
| **It holds on someone else's data** | On 135 datasets from Horlacher et al. 2023 that our panel does not contain, 108 proteins, their windows and their negatives: the directional ratio between their two negative-set constructions is **1.69** (95% CI 1.41 to 2.03) two-stage and **1.66** (1.40 to 1.98) cross-fitted, so it is comparable with both our 5.42 and our primary 4.84. Rebuilding their folds so no chromosome is split, which their release does not do, gives 1.74 and 1.68 with cross-fold near-neighbour leakage falling to exactly zero. Criteria for all four were committed before each was computed. It is disjoint in datasets and processing, **not** in proteins: 31 of 108 overlap ours, and excluding every shared one still gives 1.66. It is a pre-specified analysis of a held-out subset of an already-known benchmark, not a prospective search |
| **None of seven surveyed reports the baseline** | Of a targeted, non-systematic sample of seven methods and benchmarks, five build negatives by relocating genomic intervals, which leaves composition unconstrained, and **none** reports a composition-only AUROC. Seven hand-picked sources are not a systematic review, and the survey's selection rule is stated in `scripts/negative_set_survey.py` |

## Rebuild it from raw data

The offline check above needs nothing, and it verifies released results rather than
reconstructing them. Regenerating the result tables from raw data needs the window store and,
for the neural arms, a GPU. `run.sh all` covers the dinucleotide arm end to end; the GC and
bias-aware sweeps were run through `cloud/modal/` and their per-window scores are committed
rather than rebuilt by the default path. `docs/REPRODUCE.md` marks which is which per stage.

```bash
export GOOGLE_CLOUD_PROJECT=your-project
./run.sh preflight        # spends nothing, gates everything
./run.sh all              # pauses before every paid stage
```

Full procedure: **[docs/REPRODUCE.md](docs/REPRODUCE.md)**.
Why dataset counts differ between analyses: **[docs/PANELS.md](docs/PANELS.md)**.

Cost: **~$20 of real money** was spent on the published run and **~$60** is the forecast for a
rerun without credits. One table, with what is measured separated from what is forecast:
**[docs/COST.md](docs/COST.md)**. Every paid stage asks first.

## Layout

```
manuscript/     the paper and its figures
scripts/        one analysis per file; each writes a table under results/tables/
src/rbp/        the library the scripts import
tests/          855 tests, no network or cloud; 2 modules need torch
config/         params.yaml (the study's settings), golden.yaml (expected values)
results/tables/ every number in the paper (SCHEMA.md documents the columns)
data/evidence/  per-window out-of-fold scores for all three model classes
```

## Design rules

1. **Verification is a stage**, not an afterthought. Reproducibility that is not checked is not
   reproducibility.
2. **No hardcoded project id.** Everything resolves through `rbp.utils.cloud`, and a test fails
   the build if a literal reappears. The test watches a pattern, not one historical name: it
   watched only the old project id for a while, and the current one duly reappeared in an
   argparse default underneath it.
3. **The panel is an artefact, not a flag.** Written once, and committed three ways so a
   reader never has to take it on trust: `results/tables/supplementary_table_s1.csv` is the
   study panel of record, with ENCODE accession, experiment and an `in_three_arm_panel` flag
   per row; `config/panel_final_{cell}_{arm}.tsv` carries each arm's own membership and pair
   counts; and `manifest/study_panel.tsv` in the derived bucket is what the cloud stages read.
   `tests/unit/test_panel_is_committed.py` checks the first two against the committed
   per-window scores, because for a while the bias-aware arms had no panel file at all and
   their membership was whatever directories existed on one laptop.
4. **Task counts come from manifests**, never typed by hand.
5. **Completion markers are written last**, so an interrupted stage redoes its work rather than
   being skipped.

## Licence and citation

Code and derived data under MIT and CC BY 4.0 respectively; see `LICENSE`. Intermediate window
tables containing genomic sequence are not redistributed and are regenerated from the ENCODE
accessions in Supplementary Table S1.
