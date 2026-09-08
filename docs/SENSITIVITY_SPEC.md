# Pre-specification: the class-ratio sensitivity and the k-mer capacity ladder

**Written 2026-09-07. Committed before either analysis was run.**

Both analyses were proposed on 2026-09-06 in a costed plan and they were
dropped when the machine allocated for that run halved. They have since been described in the
manuscript as specified but unrun, which an audit correctly said invites the obvious reviewer
request. This document fixes the settings; the results are not known at the time of writing.

**They cost nothing.** The plan assumed the window tables would have to be pulled from a bucket.
They are on this machine: `rbp-store/processed/{gc,dinuc,neg2}` resolve through symlinks to
`rna-binding-proteins/data/processed`, `data/processed_dinucmatch` and a local copy, 457 MB, 187
plus 94 `dataset.tsv` files, each carrying `seq_rna`, `label` and `fold`. So both analyses are
laptop CPU work at zero spend, not the roughly \$1 of GCP the plan budgeted.

---

## B. Class-ratio sensitivity

**Question.** Every nested fit in this paper is estimated at roughly 1:1 positives to negatives.
A genome-wide screen faces a far smaller positive fraction. AUROC is prevalence-invariant, but
the fitted logistic is not: the intercept and, at fixed penalty, the effective shrinkage on each
coefficient both move with class balance. So the question is whether the **arm ordering** and
the **span** survive a change in class ratio, not whether the AUROC does.

**Ratios.** 1:1 (the published balance, as the control), **1:2**, **1:4** and **2:1**, expressed
as positives:negatives.

**How the ratio is reached.** The committed windows are approximately 1:1, so no ratio can be
reached by adding rows. It is reached by subsampling, **within fold**, so the fold structure and
its blocking are untouched:

* for 1:2 and 1:4, keep every negative and subsample **positives** to `round(n_neg / 2)` and
  `round(n_neg / 4)` within each fold;
* for 2:1, keep every positive and subsample **negatives** to `round(n_pos / 2)` within fold.

Subsampling is without replacement, seed **20260907**, one `default_rng` per (dataset, arm,
ratio) seeded from that constant plus a stable hash of those three, so the draw is reproducible
and does not depend on iteration order.

**Exclusion.** A dataset is excluded from a ratio if any fold would be left with fewer than 20
rows of either class, since a within-fold AUROC is meaningless there. Exclusions are counted per
ratio and reported. No dataset is excluded on the basis of its result.

**Estimand.** The same nested contribution as the paper, `AUROC(composition + score) -
AUROC(composition)`, pooled out of fold, with the 4-mer base model, in all three arms. Panel
means are means over datasets. The span is `max(panel mean) / min(panel mean)` across the three
arms, reported per ratio.

**What is reported, whatever it says.** The panel mean per arm per ratio, the span per ratio,
whether the ordering `dinucleotide > GC > bias-aware` holds at each ratio, and the exclusion
counts. **If the ordering breaks at any ratio, that is the finding and it goes in the
manuscript.** No protein-clustered interval is computed for the subsampled ratios, because the
subsample adds a variance component this design does not separate; the point estimates and the
ordering are what the analysis is for, and that limitation is stated rather than papered over.

---

## C. The k-mer capacity ladder

**Question.** The paper reports the outer-fold channel at one point, `k = 4` (and `k = 2` as the
floor model), and asserts a direction for how it scales with base-model capacity. This measures
the shape instead of asserting it.

**Ladder.** `k = 2, 3, 4, 5, 6`, all three arms, both estimators, so the channel at each rung is
`two-stage minus cross-fitted`. The cross-fitting procedure is `scripts/cross_fitting.py`
unmodified: ten complement fits per dataset per arm, one per unordered fold pair.

**Dataset sample, and the rule is fixed here rather than the count.** A `k = 6` design has 4096
columns, so a rung costs far more than `k = 4` and the full 94-dataset panel at all five rungs
does not fit in a sensible local budget. The rule: **systematic sampling by pair rank**, take
every `m`-th dataset from the panel ordered by pair count, with `m` the smallest integer for
which the projected runtime is under **90 minutes**, projected from a timed two-dataset
measurement. This is the same selection rule `scripts/select_panel.py` uses and for the same
reason: it spans the size range instead of taking the biggest or the smallest. The realised `m`
and dataset count are reported, and the sampled set spans the panel's pair range.

**Exclusion.** A rung is skipped for a dataset if a fold lacks both label classes. Counted and
reported.

**What is reported, whatever it says.** Per rung and arm: the two-stage contribution, the
cross-fitted contribution, and the channel between them; and per rung the three-arm span under
each estimator. The pre-committed reading is that the channel should grow with `k`, because a
higher-capacity base model has more scope to encode the withheld fold. **If it does not grow, or
is not monotone, that is reported as measured and the manuscript's assertion about direction is
withdrawn.**

**What this cannot do, stated before it is run.** It is a diagnostic and **not a substitute for
neural cross-fitting**, and no extrapolation from this ladder to the CNN or SpliceBERT will be
made or reported. A k-mer logistic and a fine-tuned transformer differ in more than the number
of parameters. The ladder bounds nothing for either neural model; that limitation stays exactly
as the Methods already state it.
