# Costed plan: redrawing the composition-matched negatives

**Status: A ran on 2026-09-06. B and C did not.** Approved by the author, executed on one VM,
about $2. The plan below is left as it was written, because the differences between what was
proposed and what happened are the useful part of the record. They are set out here first.

| | proposed | actual |
|---|---|---|
| machine | `c2d-standard-16` | `c2d-standard-8`; the binding quota is the GLOBAL `CPUS-ALL-REGIONS` limit of 12, not the regional 32 this document cites |
| analyses | A, B and C bundled | **A only.** B and C were dropped when the machine halved and are still unrun |
| dinucleotide draw | ~113 min | **~10 min.** The estimate came from a smoke test on two datasets at the bottom of the panel's size distribution, and work scales with pairs |
| wall clock | 6.0 h | 4 h 51 m, both arms at ten draws, self-deleted two hours inside its deadline |
| cost | $4.51 | about $2 |

**Result.** Both composition-matched arms carry ten draws. Between-draw SD 0.00041 (GC) and
0.00065 (dinucleotide) against between-protein SE 0.00324 and 0.00651, so propagating the draw
widens the panel-mean interval by 0.81% and 0.49%, against 3.6% on the bias-aware arm. All three
protocols now have a draw estimate.

**Corrected after an external audit.** The first version of this section reported 0.08% and
0.05%, because the summariser propagated SD/sqrt(10), the standard error of the mean over ten
draws. The published estimate uses ONE draw, so the omitted component is the SD across draws.
The error understated the widening tenfold and disagreed with `negative_draws.py`, which had
used the SD directly for the bias-aware arm all along. The conclusion is unchanged: the draw
remains the smaller term by a factor of 8 (GC) and 10 (dinucleotide). Tables: `results/tables/redraw_composition.csv` and the
per-draw files under `results/tables/draws/`. Gated by `verify_redraw_composition`.

**What this does NOT do**, and the reason is a design property rather than a defect: the
committed window tables hold RETAINED pairs, so a redraw starts from positives the matcher had
already filtered. It measures draw variability conditional on those positives. At the published
seed the panel means agree to 0.000569 (GC) and 0.000013 (dinucleotide) while individual
datasets differ by up to 0.0213 and 0.0309.

**Five launches, four of them operator error**, all now gated in the startup script: an
unpushed commit left the clone silently on `main`; a missing dependency found after staging
4 GB; the GC arm computing its region pools twice; `OMP_NUM_THREADS=1` making `$(nproc)` return
1 so the parallelism silently switched off; and a hardcoded `--zone` in the self-delete that a
capacity shortfall in that zone exposed.

---

This is the largest remaining gap in the paper's uncertainty statement. Methods discloses that
every interval is conditional on one negative draw per protocol, and that the redraw is the
largest of the six unpropagated components. Five redraws exist for the bias-aware arm only
(`negative_draws.csv`): between-protein SE 0.00135, between-draw SE 0.00037, and propagating the
draw widens the interval by 3.6%. The GC and dinucleotide arms have one draw each and no
estimate of that component at all, so the paper currently says "the ordering does not depend on
the draw" for one arm of three and must not say it for the others.

## What would be run

Three analyses, one VM, one dataset pull. They are bundled because they need the same inputs and
splitting them would triple the setup cost for no scientific gain.

### A. Multiple negative draws, GC and dinucleotide arms, 4-mer only

Ten independent redraws per arm. For each draw: resample the candidate pool from the genome at
the configured `pool_multiple`, rerun the matcher, refit the 19-column composition baseline and
the 4-mer, and recompute the pooled out-of-fold nested contribution on the study folds, which are
frozen and are not redrawn. Positives are unchanged; only the negatives differ.

**Why ten and not five.** The quantity being estimated is a standard deviation across draws. Its
own relative standard error is approximately $1/\sqrt{2(k-1)}$: 35% at $k=5$, 24% at $k=10$, 17%
at $k=20$. Ten is the point where the component can be reported with an interval rather than as a
point estimate, and it is twice what the bias-aware arm already has, so the composition arms will
not be the weaker evidence. Twenty would be better and costs roughly twice as much; if the run
comes in under budget the marginal draws are the first thing to spend the remainder on.

**The neural models are not part of this.** New negatives mean new rows, and new rows would need
new CNN and SpliceBERT sweeps at roughly $57 per full sweep across three arms. This is 4-mer
only, which is where the headline estimand lives.

### B. Class-ratio sensitivity

The nested fit is estimated at 1:1 in every arm, and a genome-wide screen faces a far smaller
positive fraction. AUROC is prevalence-invariant but the fitted logistic is not. Subsample rows
from the committed windows to 1:2, 1:4 and 2:1 and recompute. Cannot be done from the released
per-window score files, which carry no sequence, so it needs the same window pull as A.

### C. The k-mer capacity ladder for the estimator floor

Cross-fit at $k = 2, 3, 4, 5, 6$ in all three arms and report how the outer-fold channel scales
with base-model capacity. **This is a diagnostic sensitivity analysis and is not a substitute for
neural cross-fitting.** It cannot be: a k-mer logistic and a fine-tuned transformer differ in
more than capacity, and no extrapolation from the ladder to the CNN or SpliceBERT is proposed or
would be reported. What it can show is the direction and rough shape of the relationship, which
is currently asserted from one point.

## Machine, region and data

| item | value |
|---|---|
| region | `us-central1`, because both buckets are `US-CENTRAL1` and same-region reads are free |
| machine | `c2d-standard-16` (16 vCPU, 64 GB), on demand |
| disk | 200 GB `pd-balanced` |
| quota | `CPUS` in us-central1: 32 limit, 0 in use. 16 fits with headroom |
| **not spot** | `PREEMPTIBLE_CPUS` limit is 0 in this region. A preempted run at hour five would cost the whole run, and the saving is a couple of dollars |
| inputs | `gs://<derived>/processed/{gc,dinuc}/` (484 MB + 485 MB), `gs://<raw>/GRCh38.primary_assembly.genome.fa` (3.15 GB), `gs://<raw>/peaks/` (21 MB), and `data/interim/regions.pkl` |
| egress | **zero**, in region. Pulling the same data to the laptop instead would be about $0.12, which is affordable but puts the CPU work on a laptop, against the standing rule |

## Cost, and the hard stop

| line | rate | hours | cost |
|---|---|---|---|
| c2d-standard-16 | $0.752/h | 6.0 | $4.51 |
| 200 GB pd-balanced | $0.027/h | 6.0 | $0.16 |
| storage reads, same region | $0 | | $0.00 |
| egress | $0 | | $0.00 |
| **expected total** | | | **$4.67** |
| **hard maximum** | | 9.0 | **$7.00** |

**The hard stop is $7.00 and it is enforced by three independent mechanisms, because a budget
that depends on remembering to check it is not one.**

1. A `--max-runtime 9h` shutdown timer set on the instance at creation, so the VM deletes itself
   whether or not anything else works.
2. The existing `cloud/killswitch/`, whose budget threshold is lowered to $10 for the duration.
   It detaches billing, and `cloud/killswitch/main.py` is explicit that a firing should be
   treated as data loss until checked.
3. `cloud/cost.sh` polled every 15 minutes. It exits non-zero on a failed query rather than
   reporting zero spend, which is the bug it was written to fix.

**Stopping rules.** Stop and report without finishing if: elapsed time passes 7 hours; measured
spend passes $6; the first arm's ten draws disagree with the committed single draw by more than
0.01 in panel mean, which would mean the redraw is not reproducing the published construction and
the code is wrong rather than the estimate being informative; or any draw fails to produce 94
datasets.

## Intermediate output, failure recovery, cleanup

Each draw writes its per-dataset table to `gs://<derived>/redraw/{arm}/seed-{n}.csv` **as it
completes**, not at the end. A run killed at hour eight keeps every completed draw, and the
analysis in D below runs on whatever number of draws exist, with $k$ reported. Restart is
therefore idempotent: the driver skips a seed whose output object already exists.

Cleanup is `gcloud compute instances delete`, run explicitly and then verified by listing
instances in the region, not assumed from the command's exit code. The killswitch threshold is
restored afterwards. `cloud/cost.sh` is run once more the following day, because billing reports
lag by up to 24 hours and the figure at teardown is not the final one.

## Exact outputs and the analysis on them

- `results/tables/negative_draws_gc.csv`, `negative_draws_dinuc.csv`, and per-dataset tables:
  for each arm, the ten draws' panel means, the between-draw SD with its own interval, the
  between-protein SD for comparison, and the percentage by which propagating the draw widens the
  headline interval. Exactly the format `negative_draws.csv` already uses for the bias-aware arm,
  so the three arms become comparable.
- `results/tables/class_ratio.csv`: the nested contribution at 1:2, 1:4, 1:1 and 2:1 per arm, and
  whether the arm ordering survives each.
- `results/tables/capacity_ladder.csv`: the cross-fitted contribution and the outer-fold channel
  at $k = 2, 3, 4, 5, 6$ per arm.
- Golden entries and `verify.py` assertions for all three, gated like everything else.

**What the paper would then say, and what it would not.** Methods' list of six unpropagated
components loses its largest member for all three arms rather than one, and the sentence about
the ordering not depending on the draw becomes supportable for the composition-matched arms
instead of being confined to the bias-aware one. It would not become an unconditional interval:
the panel selection, the fold partition, the model class and the neural initialisation are still
held fixed, and a redraw of the negatives does not touch any of them.

## Recommendation

**Worth running, and the decision is the author's.** It closes the paper's own largest disclosed
gap in uncertainty, it costs under $7 against a $10 ceiling with three independent stops, and
nothing about it can change a published number: the redraws are a variance component reported
beside the existing estimate, not a replacement for it. The risk is that it produces a wider
interval, which is the point of running it.

The thing it cannot fix, and which no affordable experiment can, is neural cross-fitting at about
$570. That stays undone and disclosed.
