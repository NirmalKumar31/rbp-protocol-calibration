# Correction 1 to the external sensitivity analysis

**Written 2026-09-07, after the results it corrects were computed and reported. This is a
post-hoc correction and is labelled as one.**

`docs/EXTERNAL_BENCHMARK_AMENDMENT.md` was written before its analyses ran and its outcomes were
unknown at that commit. **This document is not that.** The defect below was found by an external
audit reading the committed code after the results were in the manuscript. Nothing here was
predicted; the correction is being made because the implementation did not do what the
specification and the manuscript said it did.

The original values are preserved in this file and in git history. They are not quietly
replaced.

---

## Defect 1: the two arms were given different chromosome-to-fold maps

**The specification** (`EXTERNAL_BENCHMARK_AMENDMENT.md` section 3) says positives and negatives
are refolded together so that a chromosome appears in exactly one fold, and that both arms are
refolded "by the same rule, so the two arms remain comparable". The manuscript says the fold
partition is held fixed while the negative construction varies. That is the entire logic of the
comparison.

**What the code did.** `scripts/external_sensitivity.py` called `chrom_folds(d.chrom.values)`
once per arm, inside the loop over `negative-1` and `negative-2`. `chrom_folds` balances fold
sizes by the window counts of whatever it is given. The two arms share their positives but have
different negatives, so the counts differ, so the greedy assignment diverges.

**Measured, not assumed.** Reconstructing both maps for all 135 datasets:

* identical arm maps: **0 of 135**
* different arm maps: **135 of 135**
* examples: 20 of 24 shared chromosomes remapped for `AARS:K562`, 16 of 24 for `ABCF1:K562`,
  24 of 25 for `AGGF1:HepG2`

Each arm was individually chromosome-blocked, which is why the verifier passed: it asserted that
every chromosome falls in one fold, which was true of each arm separately. It never asserted the
two maps were the same. **The gate checked a weaker property than the prose claimed**, which is
the recurring failure mode in this repository and is recorded as such.

**Consequence.** The published chromosome-blocked figures, 1.743 two-stage and 1.682
cross-fitted, mix a fold-design difference into a comparison that was supposed to isolate
negative construction. They are withdrawn as the same-fold sensitivity.

**The repair.** One map per dataset, built from an **arm-invariant source**. The positives are
identical in both arms, verified rather than assumed, so the map is built from the positive
chromosome counts alone. Chromosomes appearing only among negatives are rare (zero or one per
dataset) and are assigned afterwards in name order to the currently emptiest fold, taking the
union over both arms so the result cannot depend on which arm is processed first. The identical
map is then applied to both arms and the identity is asserted, not hoped for.

## Defect 2: the "same-strand" leakage metric discarded strand

**What the code did.** `scripts/horlacher_arm.py`'s `windows()` reads the BED6 strand field and
uses it to reverse-complement the sequence, but does not return it. `external_sensitivity.py`
then took its fallback branch, `np.full(len(d), "+")`, and computed a neighbour metric with
every window on the plus strand. That value, 0.003558, was published and described as a
same-strand rate. It is not one.

**Consequence.** The direction of the conclusion is unaffected and if anything strengthened,
because a genuine same-strand rate is smaller than a strand-agnostic one. But the number, the
population label and the denominator wording were all wrong together.

**The repair.** `windows()` preserves strand. The silent `"+"` fallback is removed and a missing
strand column is a hard failure, because a same-strand statistic computed without strand is not
a weaker measurement, it is a different one. The metric is computed for both arms rather than
`negative-1` alone, and the denominator is named exactly: windows having an eligible same-strand
neighbour within 1 kb.

## Defect 3: the no-shared-protein sensitivity was never propagated

It existed only for the original supplied-fold two-stage analysis. It is now computed for every
estimator and fold combination and gated.

---

## The original values, preserved

| cell | R as first published | 95% CI |
|---|---:|---|
| supplied folds, two-stage | 1.6912 | 1.4121 to 2.0264 |
| supplied folds, cross-fitted | 1.6640 | 1.3993 to 1.9837 |
| chromosome-blocked, two-stage | **1.7428** | 1.4412 to 2.1218 |
| chromosome-blocked, cross-fitted | **1.6822** | 1.4087 to 2.0200 |
| cross-fold near-neighbour, supplied folds | **0.003558** | strand-agnostic, mislabelled |

The two supplied-fold rows are unaffected by either defect: they use the deposit's own folds and
do not call `chrom_folds` at all. The two bolded chromosome rows and the leakage value are the
ones this correction regenerates.

**What is expected to change and what is not.** The auditor's own read-only diagnostic, using a
single shared map, obtained 1.7281 two-stage and 1.6694 cross-fitted, both still above the
protocol thresholds. That figure is a target to sanity-check against and **not** a value to
reproduce by construction: their map was built from the shared positives plus the pooled arm
windows, and this repair builds it from the positives alone, so exact agreement is not expected
and would be suspicious. Whatever the repaired analysis returns is what is reported, including
if a cell drops below the support thresholds.
