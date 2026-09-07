# Amendment 1 to the external-benchmark protocol

**Written 2026-09-07. Committed before either analysis below was run.**

This amends `EXTERNAL_BENCHMARK_PROTOCOL.md` (committed `e76a80c`, 2026-09-06 13:20:29 -0600).
**That document is not edited.** Its only value is as a record of what was fixed at `e76a80c`,
and a record corrected in place is not a record. Read the two together.

The outcomes of the two analyses specified here are **not known at the time of writing**. What
is known is the two-stage, supplied-fold result already published: span 1.69-fold (95% CI 1.41
to 2.03) over 135 datasets. Nothing below has been computed.

---

## 1. Why this amendment exists

An external audit identified three defects in the original protocol as applied. All three are
real. None is repaired by rewording the original document.

**D1. The estimand's failure rule was ill-posed.** The protocol's estimand is
`max(mean arm) / min(mean arm)`, which is bounded below by 1 by construction. Its stated
falsification rule is "fails to replicate if the confidence interval contains 1.0". On a
quantity that cannot go below 1, that interval essentially always contains 1 at its lower
limit, so the rule was close to unsatisfiable. `scripts/external_replication.py` already
computes the verdict on a **pre-labelled** `negative-1 / negative-2` ratio instead, which can
fall below 1 if the arms swap. That substitution was made after the data were inspected. It did
not change the reported number, because in 4000 draws the ordering never swaps, but the formal
rule was changed post hoc and must be re-fixed prospectively. Section 2.

**D2. The fold criterion was not checked before acceptance, and then was satisfied by a
different instrument.** Eligibility criterion 4 requires "a defined train/test partition, or
coordinates sufficient to build a chromosome-blocked one", and the falsification rule says to
report **indeterminate** if a chromosome-blocked partition cannot be made. Horlacher's supplied
partition is not chromosome-blocked: every chromosome appears in every fold. It was accepted
"in substance" on a post-hoc near-neighbour diagnostic (199 of 1,323,516 same-strand pairs
within 1 kb cross a fold, 1.5e-4), which is a finer instrument and a genuinely reassuring
measurement, but it is not the criterion that was written. The literal criterion is testable
directly, because the deposit ships coordinates. Section 3.

**D3. The external analysis uses the two-stage estimator while the internal primary estimator is
cross-fitted.** So the external number is comparable with the internal 5.42 and not with the
internal primary 4.84. The manuscript says so, but a like-for-like external figure is
computable from the same inputs. Section 4.

D2 and D3 are addressed by running the analyses. D1 is addressed by fixing the estimand and rule
in this document, before those analyses are run.

---

## 2. The directional estimand and the decision rule, stated mathematically

Let the two released negative constructions be indexed `a = 1` for `negative-1` and `a = 2` for
`negative-2`. For dataset `d` in the analysed set `D`, let

    G(d, a)  =  AUROC( composition features + model score )  -  AUROC( composition features )

pooled out of fold, exactly as in Methods, with the composition block and the model class fixed
across `a`. Define the panel means and the **directional** estimand

    M(a)  =  (1/|D|) * sum over d in D of G(d, a)

    R     =  M(1) / M(2)

`R` is directional: the numerator and the denominator are fixed by the label in the deposit,
before any value is computed, so `R < 1` is attainable and means the ordering is the reverse of
ours. This replaces `max/min`, which cannot express that outcome. The descriptive
`max(M)/min(M)` continues to be reported for continuity with the published figure, and is no
longer the quantity any decision is taken on.

`M(2) <= 0` makes `R` undefined; that case is reported as **indeterminate**, not as a number.

**Interval.** 95% percentile bootstrap over 4000 draws, **resampling proteins** (not datasets),
recomputing `M(1)`, `M(2)` and `R` within each draw. Protein is the cluster because a protein
can contribute two datasets, one per cell line. Seed **20260907**, fixed here.

**Decision rule, fixed now, for each analysis in sections 3 and 4 separately:**

| outcome | condition |
|---|---|
| **supports** | `R` point estimate `>= 1.5` **and** the interval's lower bound `>= 1.2` |
| **fails** | the interval contains 1.0, **or** the point estimate is `< 1.0` |
| **weak** | neither of the above, i.e. directionally consistent but under the thresholds |
| **indeterminate** | `M(2) <= 0`, or fewer than 20 datasets survive |

The 1.5 and 1.2 thresholds are carried over unchanged from the original protocol, which fixed
them before the 135-dataset complement was scored. Only the quantity they apply to is changed,
from a bounded-below-by-one maximum ratio to a directional one. On the published sample this is
not a loosening: the published `max/min` span is 1.69 with lower bound 1.41, and `R` equals it
whenever `M(1) > M(2)`.

**"Fails" is a real possibility and will be reported as such.** If either analysis in sections 3
or 4 fails or comes out weak, that outcome goes into the manuscript with the same prominence as
a supporting one, and the claim is downgraded accordingly. This is written down before the
result is known so that it cannot be renegotiated after.

---

## 3. Chromosome-blocked refold, specified before running

**Construction.** For each of the 135 analysed datasets, discard the deposit's supplied fold
labels and assign folds by chromosome:

1. Take the chromosome of each window from the deposit's BED coordinates.
2. Order the distinct chromosomes present in that dataset by descending window count, breaking
   ties by chromosome name ascending, so the assignment is a deterministic function of the data
   and not of dictionary order.
3. Greedily place each chromosome into whichever of the 5 folds currently holds the fewest
   windows. This is the standard largest-first bin-packing heuristic and it balances fold sizes
   without ever splitting a chromosome.
4. A dataset is **excluded** from this analysis if it has fewer than 5 distinct chromosomes, or
   if any fold ends up with only one label class, since neither can produce an out-of-fold
   AUROC. Exclusions are counted and reported; they are not silently dropped.

**Positives and negatives are refolded together**, by the window's own coordinate, so that a
chromosome appears in exactly one fold across both classes. Both `negative-1` and `negative-2`
are refolded by the same rule, so the two arms remain comparable.

**Verification, not assertion.** The analysis emits, as gated rows: the number of datasets
whose chromosomes each fall in exactly one fold (must equal the analysed count), and the
cross-fold same-strand-within-1-kb neighbour fraction under the new partition, which must be
**no larger** than the 1.5e-4 measured on the supplied folds. If the refold does not reduce or
match that leakage measure, the refold is worse than what it replaces and that is reported.

**Both partitions are reported side by side.** The supplied-fold result is not withdrawn; it is
what the deposit's own authors intend and it remains the primary external figure. The
chromosome-blocked figure is the sensitivity that answers the literal protocol criterion.

---

## 4. External cross-fitting, specified before running

**Procedure**, identical to `scripts/cross_fitting.py` and transported without modification.
For a dataset with 5 folds and design matrix `X` of 4-mer counts:

1. `published[r]`: the score for row `r` from the base model trained on every fold except `r`'s
   own. This is what the two-stage estimator uses and what the published external number used.
2. For each outer test fold `i`, the covariate for a row in any other fold `j` comes from the
   base model trained on the complement of `{i, j}`. For a row in fold `i` itself, `published`
   is already clean and is used. The complement models are symmetric in `i` and `j`, so there
   are 10 of them per dataset per arm, not 20.
3. `G_cf(d, a)` is then computed exactly as `G(d, a)`, with the cross-fitted covariate in place
   of the published one, and `R_cf = M_cf(1) / M_cf(2)`.

**Cost.** A 4-mer logistic fit on one of these datasets is milliseconds. A 2-dataset smoke test
of the existing two-stage path took 4.5 seconds wall clock including genome load. Cross-fitting
multiplies the fit count by roughly 3 and adds 10 complement fits per dataset per arm, so the
whole 135-dataset run is expected to be **under an hour on one laptop core, at zero cost**. No
GPU, no cloud, no spend. If it exceeds two hours it will be stopped and reported as blocked
rather than moved to paid compute without approval.

**Both estimators are reported.** The two-stage external figure is retained as the comparability
analysis, because it is what the surveyed literature computes and what the published number is.
The cross-fitted figure is the like-for-like comparison with the internal primary 4.84.

**Expected direction, recorded so it can be wrong.** Internally, cross-fitting moved the span
from 5.42 to 4.84, a 12% reduction, because it closes an outer-fold channel that inflates the
bias-aware arm most. If the same mechanism operates externally, `R_cf` should be somewhat
smaller than 1.69 and could fall below the 1.5 threshold, in which case the cross-fitted
external analysis is **weak** rather than supporting, and that is what will be reported. It is
also possible that it moves the other way. Writing the expectation down is the only way the
outcome can contradict it.

---

## 5. Common specification for both analyses

**Analysed set.** The same 135 datasets: Horlacher ENCODE datasets whose `protein:cell` key is
absent from our 95-dataset panel. That set is a deterministic function of the deposit and our
committed panel, and it is unchanged from the published analysis.

**Missing data.** A window whose sequence cannot be recovered from GRCh38 is dropped, and the
count of dropped windows is reported per analysis. A dataset that loses its second label class
is excluded and counted. No dataset is excluded on the basis of its result.

**Minimum sample.** Fewer than 20 surviving datasets in an analysis makes that analysis
**indeterminate**, as in the original protocol.

**Seeds.** Bootstrap seed 20260907 for both. The fold construction in section 3 is deterministic
and uses no seed.

**Genome.** GRCh38 primary assembly, the same file listed in `results/tables/raw_inputs.csv`
with its recorded MD5.

**What gets committed.** Every row produced, whatever it says, into
`results/tables/external_replication.csv` and its per-dataset companion, with golden entries and
`verify.py` assertions, gated like everything else. Tests must fail if the fold criterion, the
estimator, the thresholds or the dataset count is violated.

**What this amendment cannot buy.** It does not make the original search prospective. The
benchmark was known and its 45-dataset intersection with our panel had been scored before the
protocol was written. These two analyses are **pre-specified sensitivity analyses on an
already-analysed external dataset**, and the manuscript will describe them as exactly that. The
strongest honest label for the whole external component remains a pre-specified analysis of a
held-out subset of an already-known external construction.
