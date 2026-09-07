# Prospective protocol: searching for a second external benchmark

**Written and committed before any candidate benchmark was searched for, opened, downloaded or
scored.** The commit that adds this file contains no result. That ordering is the only thing
that makes what follows worth more than a post-hoc rationalisation, and it is also the only
evidence of it: this is a repository commit, not a registry entry, and it establishes the order
of two events and nothing stronger. See the corresponding paragraph in Methods.

## Why a second benchmark, and which claim it is for

The paper makes two claims that are routinely conflated and must not be here.

**Claim A, the central one.** A model's measured contribution over a composition baseline
depends strongly on how the negative set was constructed, so contributions measured under
different protocols are not comparable without the corresponding baseline. This is established
within our panel across three protocols and three model classes.

**Claim B, the weaker one.** Harder apparent discrimination produces larger measured
contribution: the directional relation between a protocol's composition baseline and the
contribution it leaves.

**Claim B does not replicate on the one external benchmark we already have**
([Horlacher et al. 2023](https://doi.org/10.5281/zenodo.10600977)): the pooled relation there is
null, and the paper says so. This search is for a benchmark that can test **Claim A**, which the
Horlacher data cannot, because testing Claim A needs a dataset released under **more than one
negative-set construction for the same positives** and Horlacher releases one.

We are explicitly not searching for a dataset that would rescue Claim B. If one is found that
bears on Claim B, it is reported, but Claim B stays labelled as not replicating unless a
qualifying dataset supports it under the criteria below.

## Eligibility criteria

A candidate qualifies only if **all** of the following hold. These are conjunctive and are
checked in this order, so that the cheapest disqualification is found first.

1. **Independent of this study.** Not derived from ENCODE eCLIP K562/HepG2 peak calls that
   overlap our 94 datasets. A benchmark built from the same peaks under different processing is
   not independent and is excluded.
2. **At least two negative-set constructions over the same positive set**, released as data
   rather than described in prose. This is the criterion that makes the benchmark able to test
   Claim A at all, and we expect it to be the one that fails.
3. **Sequences are recoverable** for both positives and negatives, either released directly or
   reconstructible from coordinates plus a named genome assembly. Without sequence there is no
   composition baseline and therefore no contribution to measure.
4. **A defined train/test partition, or coordinates sufficient to build a chromosome-blocked
   one.** A random split is disqualifying, because a leaky partition inflates every contribution
   and would confound the comparison with our chromosome-blocked panel.
5. **At least 20 RBP datasets**, so that a panel mean has a resampling distribution and a
   clustered interval is meaningful. Fewer than 20 is reported as found-but-underpowered rather
   than as a test.
6. **Public and redistributable**, with a licence permitting reuse, and reachable without an
   access request.

## Exclusions

Excluded regardless of the above: in-vitro-only assays (RNAcompete, RNA Bind-n-Seq) where the
negative is a synthetic library rather than a genomic window, because the composition baseline
means something different there; benchmarks whose only negatives are dinucleotide shuffles of
the positives, since our own shuffled arm already shows that construction pins the baseline at
0.5 and it is a degenerate case rather than a second protocol; and any dataset for which the
negative construction cannot be determined from the release.

## Primary hypothesis and estimand

**Hypothesis.** On a qualifying external benchmark, the nested contribution of a fixed model
class over a fixed composition baseline differs across that benchmark's negative-set
constructions by a factor materially greater than one.

**Estimand.** Exactly the one defined in Methods, transported unchanged:
$S = \max_a \bar{\Delta}_a / \min_a \bar{\Delta}_a$, where $\bar{\Delta}_a$ is the unweighted
mean over the benchmark's datasets of the pooled out-of-fold nested contribution in
construction $a$. The model class is the 4-mer logistic regression and the baseline is the same
19-column composition block, both unchanged from this paper, because changing the model at the
same time as the data would make a difference uninterpretable.

**Uncertainty.** Percentile bootstrap, 4000 draws, clustered on protein if the benchmark
identifies proteins and on dataset otherwise, matching this paper's procedure. The interval is
conditional on the benchmark's own fold partition and negative draws in exactly the way ours is
conditional on ours.

## Falsification criteria, fixed now

Stated as decisions, not as directions to interpret afterwards.

- **Claim A is supported** on the external benchmark if $S > 1.5$ and the lower bound of its
  95% interval exceeds $1.2$.
- **Claim A fails to replicate** if the 95% interval for $S$ contains $1.0$.
- **The result is indeterminate**, and reported as such rather than as either of the above, if
  $S$ falls between those two regions, or if fewer than 20 datasets qualify, or if the
  benchmark's fold partition cannot be made chromosome-blocked.
- The **arm ordering** is a secondary outcome and is not a falsification criterion, because the
  external benchmark's constructions need not correspond to ours and there is no reason its
  ordering should match.

No outcome of this search changes any claim already made about our own panel. A failure to
replicate on an external benchmark would restrict the scope of Claim A to eCLIP-derived panels,
and would be reported in Limitations in those words.

## What is reported if nothing qualifies

The search itself, in full: the sources queried, the candidates considered, and the criterion
each failed on. **A null search result is the expected outcome and is a finding.** If no public
RBP benchmark releases more than one negative-set construction over the same positives, that is
itself the strongest available evidence that the comparison this paper warns about cannot
currently be checked by a reader, which is a reason for the recommendation rather than a gap in
it. It will not be presented as if a qualifying benchmark were found and then dismissed.

## Search plan

Sources, in the order to be queried: the negative-set constructions cited in this paper's own
targeted survey and their data releases; Zenodo, figshare and GEO for RBP binding benchmarks;
the data-availability statements of the RBP prediction methods published since 2020 that the
survey names; and any benchmark suite advertising multiple negative constructions. Everything
found is recorded in `results/tables/external_search.csv` with its disqualifying criterion,
including the candidates that fail on criterion 1 or 2 immediately.


---

# Addendum, written after the search. The protocol above is unedited.

Recorded here rather than folded into the text, because a prospective protocol whose errors are
corrected in place afterwards is not a prospective protocol.

**The title and opening claim of this document are too strong, and an external audit was right
to say so.** "Before any candidate benchmark was searched for, opened, downloaded or scored" is
false of the benchmark this search actually landed on. The Horlacher deposit had been downloaded
and MD5-verified on 2026-08-31 (`0a79763`), its 45-dataset intersection with our panel had been
scored the same day (`3f96e8e`, `2b5843a`), and it was analysed further on 2026-09-04
(`28a5a1f`), all before this file was written on 2026-09-06.

What is true, and is what the manuscript now says, is narrower: **the criteria below were
committed before the 135-dataset disjoint complement was scored.** That complement had never
been touched, so the rule applied to it was genuinely fixed in advance. Call this a
pre-specified analysis of a held-out subset of an already-known external construction. It is
weaker than a prospective benchmark search and it is not nothing.

The opening paragraph is left standing rather than rewritten, for the same reason the rest of
the protocol is: this document's only value is as a record of what was committed at `e76a80c`,
and a record that gets corrected in place is not a record. **Read it with this addendum.**

**The fold criterion above was never checked until an audit asked, and the answer is two-sided.**
The protocol calls the result indeterminate if the benchmark's fold partition "cannot be made
chromosome-blocked". Horlacher's partition is **not** chromosome-blocked: every chromosome
appears in every fold, in all 135 datasets of the complement.

Measured rather than left there, because chromosome grouping is a coarse instrument for one
specific channel and the channel is what matters. The direct metric, the one
`fold_integrity.py` applies to our own arms, is the fraction of positives whose same-strand
neighbour within 1 kb sits in a different fold. On their coordinates:

    positives with a neighbour within 1 kb      1,323,516
    of those, neighbour in a different fold           199   = 1.5e-4
    worst single dataset                                     0.905%

So their scheme is **locus-blocked rather than chromosome-blocked**, and closes the leakage
channel more finely than ours does. The criterion is satisfied in substance and not in form, and
both halves are now gated so neither can be quoted without the other.

**The estimator transported here is the two-stage one, not the cross-fitted primary.** The
protocol says the estimand is "exactly the one defined in Methods, transported unchanged", and
at the time it was written the paper's headline estimator WAS the two-stage one. The paper has
since made the cross-fitted estimator primary, which makes the external span comparable with the
5.42 rather than with the 4.84. The abstract said "the same estimator" and now names it.
Cross-fitting the external arms is affordable and has not been run.

**A premise above is also wrong.** The protocol states that the Horlacher benchmark "releases one"
negative-set construction and therefore cannot test Claim A. It releases **two**, per fold, over
the same positives: `negative-1.fold-N.bed` and `negative-2.fold-N.bed` beside
`positive.fold-N.bed`, for all 223 ENCODE datasets in the deposit. The error was in this
repository's own reading of a benchmark it already had on disk, and it meant the paper had been
treating an available external test of its central claim as unavailable.

**Nothing else in the protocol changes.** The eligibility criteria, the estimand, the
uncertainty procedure and the falsification thresholds are as committed at `e76a80c` and were
applied without modification. In particular criterion 1 was applied strictly, which is what
sends the analysis to the 135 datasets our panel does not contain rather than to the 45-dataset
intersection that `horlacher_arm.py` already publishes.

**The search itself is in `results/tables/external_search.csv`**, with every candidate and the
criterion each failed on. The strongest remaining candidate is the Mukherjee PAR-CLIP subset of
the same deposit, which would widen scope beyond eCLIP; it is not analysed, because a different
crosslinking chemistry has different positional bias and the protocol fixes the estimand as
transported unchanged. It is recorded as the next experiment rather than as a result.

**The result is in `results/tables/external_replication.csv`** and was produced by applying the
thresholds above mechanically in code rather than in prose.

---

**SUPERSEDED IN THREE RESPECTS BY `EXTERNAL_BENCHMARK_AMENDMENT.md`, 2026-09-07.** This document
is still not edited, for the reason given at the top of this addendum. But three of its
provisions were found defective by a later audit and the amendment replaces them. Read the
amendment alongside this file:

1. **The estimand's failure rule was ill-posed.** `max/min` is bounded below by one, so "fails
   if the interval contains 1.0" was close to unsatisfiable. The amendment fixes a directional
   `negative-1 / negative-2` ratio instead, which can fall below one.
2. **Eligibility criterion 4 was never checked before the benchmark was accepted**, and was
   then satisfied by a post-hoc near-neighbour diagnostic rather than by the chromosome-blocked
   partition the criterion names. The amendment builds that partition. It succeeds on all 135
   datasets and the result holds.
3. **The estimator here is two-stage**, while the paper's primary estimator is cross-fitted.
   The amendment runs the cross-fitted external analysis, so the comparison is like-for-like.

The amendment was committed at `02a2bac` before either analysis it specifies was run. Its
results are in `results/tables/external_sensitivity.csv` and all four cells meet the thresholds
this document fixed.
