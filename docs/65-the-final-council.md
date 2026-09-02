# 65. The final council: four referees, two rounds, and what it cost

Convened 2026-08-30/31, after R1g. Four independent reviewers, each required to RUN CODE rather
than read: a statistician, an ML-methodology referee, a genomics domain referee, and a
reproducibility auditor whose only job was to forge results past the verifier.

**Every one of them found something real.** Three findings were mine, made in this session, and
one of those I had already "independently verified" using the very assumption that was wrong.

---

## What was withdrawn

| claim | status | who killed it |
|---|---|---|
| "the contrast GROWS with model capacity" | **withdrawn** | genomics (ratio scale) + ML (apparent-drop ladder) |
| "the protocol effect rises with capacity" | **withdrawn** | statistician (baseline non-invariance) |
| "the k-mer's protocol effect survives" | **withdrawn** | statistician round 2 (slope-source grid) |
| "protocol effect is +0.0188 to +0.0313" | **too narrow**, superseded | statistician (link family) |
| "identical seed, so hardware changes speed not the model" | **false** | ML (seed ordering bug) |
| strand artifact "−0.0036, 90.6% surviving" | **stale, retracted values** | genomics |

## What survived every attack, in both rounds

**The raw contrast: +0.0398 (k-mer), +0.0530 (CNN), +0.0864 (SpliceBERT).** It is a difference
of two AUROC differences on the same rows: no transplant, no link, no transportability
assumption. Both the statistician and the ML referee explicitly said to keep it as the headline.
It survived precision weighting, random-effects meta-analysis, protein clustering, trimming,
leave-worst-out, shared-positives restriction, and an expression-confound control.

Also clean, checked and not broken: fold discipline (0 leakage in 3.59M rows across 187
datasets per arm), the training code being byte-identical between arms, the binormal
approximation, and the summary tables being genuine arithmetic on their evidence.

---

## THE BLOCKER: the R1g evidence table is forgeable end to end

The reproducibility auditor replaced every per-dataset gain with a closed-form analytic
sequence bearing no relation to any protein, cell line or GPU run, keeping only the real
`comp_*`, `*_full_*`, `coverage_*`, `n_*` columns, regenerated the summary with `--from-cache`,
and got **314/314**. Forged SpliceBERT contrasts read as a visible monotone ramp
(0.09897, 0.09963, 0.10029, ...) and nothing objected.

**Why the arithmetic gate did not stop it.** `max_summary_arithmetic_diff` verifies the summary
is arithmetic on the evidence. That is equally satisfied by evidence reverse-engineered from the
summary. It is a one-directional consistency check, and `deep_contrast_per_dataset.csv` is a
terminal artifact: hand-committed, not regenerable from anything in the repo (neither arm's
input windows are committed), with no assertion tying it to the panel, the hardware, or the 940
score files.

**Five more that passed silently:**

1. **Deleting the coverage columns** makes the gate assert 1.0. `summarise()` reads
   `d[col].min() if col in d else 1.0`. A wrong value is caught; an absent column is not.
   Same shape as past failure (b), one level down.
2. **Dataset identity is entirely ungated.** Renaming all 94 datasets to `FAKE0:XXXX`, or
   permuting the label column so every row attaches to the wrong protein: 314/314.
3. **Zeroing the `*_se_*` columns and setting `n_gc`/`n_dn` to 1**: 314/314.
4. **The weak rung can be smoothed over.** `min_step_datasets_cnn_minus_kmer: 55` is one-sided,
   so rewriting 58/94 to 94/94 passes. `golden.yaml` claims this number "is stated, not smoothed
   over". Every `positive_datasets` floor has the same defect.
5. **Every interval can be collapsed to zero width.** `ci_low = ci_high = value` passes; no gate
   bounds CI width or asserts `ci_low < value < ci_high`.

**The three fixes, in priority order.**
1. Gate `deep_contrast_per_dataset.csv` against a recorded content hash of the 940
   `scores.tsv.gz` + `metrics.json` it derives from, and cross-check its `dataset` column
   against `rehearsal_binding_gc.csv`.
2. Make `min_row_coverage_*` FAIL when the column is absent instead of defaulting to 1.0.
3. Make the one-sided `at_least` gates two-sided: inflation currently passes everywhere.

---

## The reframing, corrected: the multiplier is a property of the PROTEIN

I proposed "the protocol multiplies a model's contribution by ~2.4x-3.5x, invariant across 79
RBPs, two cell lines and three model classes". The genomics referee ran the variance
decomposition over the 259 (dataset x model) cells and **one word was wrong**:

**RECOMPUTED AND CORRECTED 2026-09-01** (`scripts/multiplier_variance.py`,
`multiplier_variance.csv`). The referee's decomposition was committed to no script and no
table -- the paper's most quotable line about R1g was sourced to a conversation. Rerun on the
262 (dataset x model) cells with both arms positive, it approximately reproduces, and **a
permutation null changes what it means**:

| factor | levels | share | null, relabelled | **excess** | p |
|---|---|---|---|---|---|
| **protein** | 79 | **64.8%** | **29.7%** | **+35.1** | < 0.0005 |
| model class | 3 | 2.8% | 0.8% | +2.0 | **0.023** |
| cell line | 2 | 0.2% | 0.4% | −0.2 | 0.49 |
| residual | | 32.8% | | | |

Panel multiplier **2.99x** (exp of the mean log; a ratio of means is a different number).

**Two things the original framing got wrong.** First, "68.9% against 1.5%" is not a comparison:
a factor with 79 levels absorbs 29.7% of the variance of *relabelled* data, so most of the gap
is degrees of freedom. The comparable statistic is each factor's excess over its **own** null.
Second, **model class is not null** (p = 0.023), which is what R1g's own withdrawal already
said when it recorded SpliceBERT's significantly lower multiplier (−0.258 [−0.369, −0.139]).
The correct claim:

> The protocol multiplier is mostly a property of the **protein**: RBP identity explains 64.8%
> of its log variance, 35 points more than a 79-level factor absorbs by chance. Cell line is
> indistinguishable from noise. Model class is small but real, and SpliceBERT's is genuinely
> lower.

Two wording cautions they gave, both worth keeping: do not launder a detected decline into
"invariance" (SpliceBERT is significantly lower, −0.258 [−0.369, −0.139]), and keep the absolute
scale alongside the ratio because they are true statements about different quantities -- what a
benchmark reader sees versus what a model developer sees.

**And the reframing repairs R1d's gated limitation.** On the 13 proteins in both cell lines,
partialling out per-protein mean total gain: the absolute contrast's cross-line correlation
collapses (+0.900 -> +0.230, p = 0.449) but the **log multiplier survives (+0.680 -> +0.580,
p = 0.038)**. Spearman with total gain across all 94: contrast +0.833, multiplier −0.221. The
absolute contrast is a rescaling of signal strength; the multiplier is not. n = 13 and the
interval touches zero, so suggestive rather than established -- but strictly better than the
limitation it replaces.

---

## The cheapest high-value experiment: a third arm, buildable from data already on disk

Horlacher's negative-2 recipe, quoted from the paper: sample negatives from *"binding sites of
other RBPs experimentally assessed in a given dataset... we only sample positives of other RBPs
which do not overlap with positives of the target RBP"*, 1:1, 101-nt windows -- identical to our
`windows.size`.

**Every ingredient is already in `rbp-store`.** For target protein P in cell C, draw 1:1
negatives from the pooled `label==1` rows of the other panel proteins in C (99 in K562, 88 in
HepG2), excluding any within `min_peak_distance` of P's own peaks. No download, no GPU, no
pipeline -- a join over files we have.

Such negatives are **100% expressed**, **100% CLIP-accessible**, **strand-correct by
construction** (each inherits the strand of the peak it came from), and **free of our composition
matcher's hyperparameters entirely**. It removes four separate objections at once and turns the
paper from a two-point contrast into a three-protocol correction factor, one of which is the
field's own, with Horlacher's published −0.065 to −0.085 AUROC drop as external calibration
against our −0.1095.

Cost: construction ~0, k-mer refit ~40 min local, optional SpliceBERT ladder ~$16 within budget.

### A disagreement between two agents, recorded

The genomics referee concluded Horlacher's negative sets are **not** downloadable, citing the
paper's data-availability statement and the GitHub repo. A dedicated availability check found
they **are**: Zenodo `10.5281/zenodo.10600977`, `samples.tar.gz`, 379 MB, CC-BY-4.0, 302
experiments x 5 folds of `positive` / `negative-1` / `negative-2` BED. It is not cited in the
paper and is discoverable only through a comment on GitHub issue #2, which is why the referee
missed it. **The referee is wrong on availability and right on the pipeline** -- their warning
not to re-run it from scratch stands (private HPC symlink, no accession list, unseeded
`sort --random-sort` fold splits). See docs/64 section 7.

Both routes are open. The local negative-2 arm is cheaper and needs nothing; the Zenodo route
gives 302 experiments and the field's own splits.

---

## Claims that could not be substantiated

- **"470/470 runs, 0 failures"**: one metrics.json predates the sweep (a local smoke test the
  resume rule adopted), so Modal ran 469; there is no GC-sweep log, and `retries=2` masks
  retried tasks.
- **"~$16.02"**: no bill or invoice artifact exists; the repo's own estimator says $30.89. The
  figure came from summing training seconds x a rate derived from a prior arm.
- **"the only difference is how the negative windows were chosen"**: the positive sets differ
  too (Jaccard 0.9164-1.0000). Tested and immaterial (contrast moves +0.0398 -> +0.0401), but
  the sentence is false as written.
- **"19.7M parameters"**: the encoder is 19,718,656; the trained model is 19,784,449 including a
  262,656-parameter pooler that is instantiated and never used.
- **"3-layer CNN"**: it has 4 learnable layers.
- **CNN provenance**: the dinucleotide CNN ran on **x86 GCP Batch CPU**, not a GPU. docs/60 says
  "the dinucleotide sweep on one A10G generation", which is wrong for that rung.

**Substantiated on re-check:** 7,089 CNN parameters, GC-arm 100% row coverage (recomputed for
all 94), 11.99 and 2.27 GPU-h, and identical seed / epochs / batch / LR / weight decay / fp32 /
fold maps across arms.

---

## The audit that audits is stale and wrong where it matters

`docs/61` lesson 10 says `audit_manuscript.py` is "6.3% saturated at four decimals, ~44% at
three". It now self-reports **11.3% / 61.9%**. Worse, saturation is computed only over
**[0.5, 1.0]**, the AUROC range -- but R1g's entire headline lives in **[0, 0.1]**, where the
4-dp grid is **61.6% occupied** and the 3-dp grid is **100% occupied**. Planting ten lies in a
copy of `docs/60`: **7 of 10 missed**, including `+0.0864 -> +0.0917`. The regex needs a decimal
point, so `470/470`, `94/94`, `7,089`, `85.4%` and `p=1.5e-05` are structurally invisible.

---

## Standing verdicts

| referee | round 1 | after fixes |
|---|---|---|
| statistician | Major revision (BLOCKER on identification) | **Minor revision** |
| ML methodology | Major revision (BLOCKER on the ladder) | not re-asked |
| genomics domain | Major revision | not re-asked |
| reproducibility | — | **BLOCKER: evidence forgeable** |

**Venue, unchanged and agreed by two referees independently:** *NAR Genomics & Bioinformatics*
(Methods) or *Bioinformatics Advances*. Not *Genome Biology* / *NAR*: those want a claim about
RNA-protein recognition and this paper has none, and its negative set is a weaker version of a
benchmark that already exists.

## Queued, all free

1. The three anti-forgery fixes above.
2. `expression_control.py` (written by the genomics referee, running, ~8 h) -- gate it.
3. Extend `cluster_intervals.py` to R1, R1b and R1f; record that R1c's 40 datasets are 40
   distinct proteins so no clustering exists there.
4. Restate 80/94 and 82/94 at the measured design effect: they become ~72/94 and ~78/94.
5. Fix the provenance sentences and the stale `audit_manuscript` figures in docs/61.
6. Re-run the CNN rung both arms on one backend with the seed fix (2/945 banked, resumable).
7. Build the negative-2 third arm.

---

# Round 3, 2026-08-31: the council that changed the thesis

Three referees: a handling editor, a hostile Referee 2, and a framing/impact adviser. All three
independently reached the same conclusion, and it was not the one the paper was making.

## What they converged on

**The baseline is the finding; the protocol label is not.** The editor put it as a positive
reframe, the framing adviser as a calibration paper, Referee 2 as a rejection reason. Same
statement.

## The three findings that changed the paper

**1. My non-monotonicity claim was false** (editor). `docs/63` called it "the guard against the
obvious rewrite" and it does not hold. Ordered by MEASURED difficulty:

| | dinuc | GC | neg2 |
|---|---|---|---|
| composition alone | 0.6274 | 0.7827 | 0.8248 |
| apparent AUROC | 0.6937 | 0.8092 | **0.8370** |
| nested contribution | **0.0663** | 0.0265 | **0.0122** |

Perfectly anti-monotone, 3 for 3. neg2 is the EASIEST discrimination, not the hardest. I had
ordered protocols by "bias-awareness", which is a subjective ranking, and called the result
irregular. **"Harder negatives reveal more" is exactly true.**

**2. The protocol label adds 1% once the baseline is known** (Referee 2). R1l had argued
protocol and baseline cannot be separated because the arms barely overlap. That was too
comfortable: Referee 2 found two places where they DO overlap.

- incremental R2 of the protocol label given the baseline: **1.00%**; of the baseline given the
  label: **11.04%**
- **the natural experiment**: neg2 usually raises the baseline, but in 27 of 94 datasets it
  lowers it -- and there the result REVERSES, −0.0212 → **+0.0028**
- matched on baseline within 0.02 AUROC, dn−gc = **−0.0087 [−0.0265, +0.0122]** against a raw
  +0.0398. The contrast does not survive matching.

This became R1n and then the thesis.

**3. Most of the magnitude is where the baseline stops** (Referee 2). A bag of 4-mers has no
positional information, so the headline is order-3/4 composition beyond order-1/2. Extending
the baseline to order 3 removes 46–56% of each arm and **64% of the R1 contrast** -- but the
**fold range survives** (2.67x → 2.62x). The magnitude is an analyst's choice; the protocol
dependence is not. This became R1o.

## What earned the title

The editor and framing adviser both flagged the first objection a referee raises: is the
5.4-fold range just the bounded AUROC scale? R1m answers it with a search rather than an
argument -- eight monotone transforms, floor **2.00x [1.67, 2.46]**, interval excluding 1, and
the transform achieving the floor divides by the protocol's own baseline rather than rescaling
the model. **This analysis could not lose**: had a transform collapsed the range, the paper
would have become "here is the protocol-independent coordinate", which is better.

## The external validation, and its limit

R1p ran the whole measurement on Horlacher's published negative sets, their folds, 45 datasets.
**The range replicates (2.38x against our 2.50x) and the baseline gradient replicates in sign
(rho −0.372, p=0.012).** But the natural experiment does NOT reverse there: −0.0409 where
negative-2 raises the baseline, −0.0184 where it lowers it, both negative. A protocol-specific
residual remains that our data does not show, **so R1n's strong form is ours only** and the
claim was weakened to "most of".

## Three desk-reject risks, all closed

1. The README advertised **three retracted claims** as its headline (composition-share doubling,
   R3 locality, R4 ClinVar transfer). An editor clicking the data-availability link would have
   returned the paper. Rewritten.
2. The coverage gate still read `if col in d else 1.0`, so deleting the evidence asserted
   perfect coverage. The anti-forgery fix I had reported as closed.
3. `0.53x` and `2.9x` are geometric means of per-dataset ratios while the panel means printed
   beside them give `0.46x` and `2.50x` -- two conventions in one document.

## The forgery, and how it was closed

The reproducibility auditor replaced every per-dataset gain with an analytic sequence unrelated
to any protein or GPU run, regenerated the summary, and passed **314/314**. The arithmetic gate
is one-directional: verifying a summary is arithmetic on its evidence is equally satisfied by
evidence reverse-engineered from the summary.

Closed by **committing both arms' 940 per-window score files (16 MB)**, so every cell's raw
pooled AUROC is recomputable inside `verify.py` from `data/evidence` alone. Plus the exact
identity `gain == full − comp`, panel-name cross-check, coverage presence, and two-sided count
gates. **Six attacks run afterwards, all caught**, including a mean-preserving forgery that only
the identity can see.

## Standing verdicts after round 3

| referee | verdict | note |
|---|---|---|
| handling editor | send out; **~64% NAR GB, ~72% Bioinformatics Advances** | "solid but minor methods paper", said as the honest answer |
| Referee 2 | **major revision**; "I could not find a fatal flaw" | attacked construction directly and failed |
| framing adviser | reframe to calibration | "the same paper, a two-hour edit" |

**Not Genome Biology (~2%).** It wants a claim about RNA-protein recognition; this paper has
none and R1f concedes it.
