# 62. How the paper was finished: every bug, every fix, every decision

**What this document is.** `docs/59` records the council rounds as narrative. This one is the
engineering and scientific record underneath them: what was actually run, what broke, why it
broke, and what the fix was. It is written to be read by someone who was not there, including
the author in six months, and it assumes no memory of the earlier documents.

It covers the period from council round 7 to the finished paper: **21 commits**, 46 scripts,
42 committed tables, 9 figures, and a verifier that grew from 112 assertions to 241.

**The one-line summary.** The project began this period with four claimed results and ended
with one, defended by five controls. Almost everything removed was removed because someone ran
code against it, not because someone thought harder about it.

---

# Part 0. Where things stood, and the shape of the problem

The study benchmarks RNA-binding-protein models on 94 paired ENCODE eCLIP datasets. Each
dataset has *positive* windows (where a protein binds) and *negative* windows (where it does
not). The negatives are not random: they are **matched** to the positives on sequence
composition, so that a model cannot win by noticing that bound regions are, say, GC-rich.

There are two matching protocols in play:

- **GC-matched** negatives, which fix only the combined G+C fraction. This is the protocol in
  general use.
- **dinucleotide-matched** negatives, which fix all 16 dinucleotide frequencies.

The paper's question is what difference that choice makes to what a benchmark measures. Not to
the headline AUROC, which obviously falls when you make the task harder, but to the **nested
contribution**: how much a sequence model adds *over and above* a 19-feature composition
baseline fitted on the same data.

    nested contribution = AUROC(composition + model score) - AUROC(composition alone)

That quantity is the paper. Everything else is a control on it.

## The four things that were claimed at the start of this period, and their fates

| claim | fate |
|---|---|
| R1, the protocol effect | **survived**, and is now the paper |
| R2, a four-model ladder | cut: Horlacher et al. 2023 published the negative-set effect across 11 methods |
| R3, pooled vs paired AUROC | cut: prior art (van Klaveren 2016; Janes & Pepe 2008), and not independent of R1 |
| R4, a ClinVar variant analysis | **retracted twice**, then cut |

---

# Part 1. R1 survives the ceiling objection

## The objection

R1's headline is a *difference of differences*: the nested contribution is +0.0265 in the GC arm
and +0.0662 in the dinucleotide arm, so the contrast is +0.0397.

Here is why that is not obviously a result. AUROC is bounded at 1. The composition baseline
sits at **0.7827** in the GC arm and **0.6280** in the dinucleotide arm. A fixed amount of real
signal added on top of a high baseline buys a *smaller* AUROC increment than the same signal
added to a low one, purely because the scale compresses near its ceiling.

The measured compression factor between those two baselines is **1.51x**. So the observed
direction of the contrast is predicted by arithmetic alone, with no protocol effect whatsoever.
That is a serious objection and it cannot be answered with a paragraph.

## The wrong fix, recorded because it is the obvious one

Somers' D. It is a rescaling of AUROC:

    D = 2 * AUROC - 1

so a nested gain on the D scale is exactly twice the same gain on the AUROC scale, and the
contrast merely doubles. **A linear rescaling cannot diagnose a nonlinearity.** This was
proposed and is worthless; it is in the paper so nobody proposes it again.

## The fix that works

Put both arms on a scale that is linear in signal rather than bounded. Under a binormal model,

    d' = sqrt(2) * Phi^-1(AUROC)

and an increment in d' is an increment in signal wherever the baseline sits. Then **transplant**
the GC arm's own d' increment onto the dinucleotide arm's baseline. That predicts what the
dinucleotide arm would show if the protocol had moved the baseline and changed nothing else.
The gap between that prediction and what the dinucleotide arm actually shows is the protocol
effect with compression removed.

    compression      +0.0083 [+0.0061, +0.0109]      21% of the contrast
    protocol effect  +0.0313 [+0.0267, +0.0363]      79%, positive in 87/94

`scripts/scale_check.py`. And the same comparison computed directly on the unbounded d' scale
gives **+0.1290** [+0.1091, +0.1499], which would be zero if compression were the whole story.

## The transplant runs both ways, and reporting one direction was question-begging

This was caught in the final council, and it matters. The same logic licenses transporting the
*dinucleotide* arm's increment onto the *GC* baseline instead. That is equally defensible and it
attributes far more to compression:

| direction and link | compression | protocol effect |
|---|---|---|
| GC increment onto dinuc baseline, probit | +0.0083 | **+0.0313** |
| dinuc increment onto GC baseline, probit | +0.0182 | **+0.0215** |
| GC onto dinuc, logit link | | +0.0288 |
| dinuc onto GC, logit link | | **+0.0188** |

The choice of direction moves the estimate further than any single interval is wide. So the
paper reports the **range +0.0188 to +0.0313** and states that each row's interval is
conditional on that transport choice. What is robust, and is gated, is that **every member of
the family keeps the sign**.

**The teaching point.** When an estimate requires a modelling choice, the honest presentation is
the family of estimates the defensible choices produce, not the member you happened to compute
first. A confidence interval conditional on a choice is not a confidence interval for the thing.

## A reversal that had to be reported

The same contrast computed on a third scale, the Firth coefficient of the standardised score in
the nested fit, **changes sign**: +1.063 in the GC arm against +0.686 in the dinucleotide arm,
dinucleotide larger in only 11/94.

A result whose sign depends on the scale is not a result, unless the reversal has a diagnosis.
It does. A logistic coefficient is identified only against the latent residual scale, so
coefficients from two fits with different total signal are not comparable (Mood 2010; Karlson,
Holm & Breen 2012). The GC task carries **1.79x** the total signal.

The fingerprint is direct and is gated: across the 94 datasets the between-arm coefficient gap
tracks each task's **total** signal at Spearman **+0.520** (p = 8e-08) and tracks the
**incremental value** it is supposed to measure at **+0.065** (p = 0.53). It is measuring the
wrong thing.

## And a "second answer" that was not one

I added a second argument: subtract what pure latent rescaling would predict, and a positive
residual of +0.1156 remains. I claimed it was stronger than the fingerprint because it is
arithmetic on the same fits rather than a correlation across datasets.

It was neither stronger nor a second answer. It is **algebraically the same statistic** as the
normalised contrast already reported:

    residual_i == dfull_gc_i * (coef_dn_i/dfull_dn_i - coef_gc_i/dfull_gc_i)

verified at max|difference| **5.0e-16** with sign agreement 94/94. And the null it is measured
against, that the coefficient scales as the first power of total signal, is a **choice**:

    exponent 0.5  ->  residual -0.1721
    exponent 1.0  ->  residual +0.1156
    exponent 1.5  ->  residual +0.5278

The sign of the "answer" is a free parameter. Both numbers stay in the tables as retraction
evidence and nothing is claimed from them.

**The teaching point.** Two derivations that use the same inputs and the same algebra are one
derivation. Check whether a corroborating statistic is actually independent before calling it
corroboration.

---

# Part 2. Own-label leakage in a baseline

A trivial baseline in the ClinVar analysis scored a variant by the *prevalence of pathogenic
variants in its 1-Mb genomic block*, leaving the variant itself out:

```python
blk = d.chrom + "_" + (d.pos // BLOCK).astype(str)
gb = d.assign(_b=blk).groupby("_b").label
tot, cnt = gb.transform("sum"), gb.transform("size")
d["prev"] = ((tot - d.label) / (cnt - 1)).where(cnt > 1, np.nan)
d = d.dropna(subset=["prev"]).drop_duplicates("vid")     # <-- too late
```

The leave-one-out subtraction is correct **only if each variant appears once**. The assignment
table carries **2.40 rows per variant**, one per dataset the variant was scored in. Subtracting
one copy leaves roughly 1.4 copies of the variant's own label inside its own block.

Fix: deduplicate before the groupby, not after. Verified safe first, 0 of 27,492 variants
disagree on label or position across their copies, so which copy survives cannot matter.

    unstratified AUROC   0.8246  ->  0.8164
    within-decile range  0.7678-0.8319  ->  0.7548-0.8173

**The teaching point.** A leave-one-out statistic on a table with duplicated rows leaks. The
duplication does not have to be a bug for the leak to be one.

---

# Part 3. The verifier, and four ways it certified nothing

The project's headline engineering artifact is `scripts/verify.py` plus `config/golden.yaml`: a
set of named numeric assertions run against committed tables. During this period it was attacked
repeatedly and failed four times.

## 3.1 Gates that vanish instead of failing

Most gates were written:

```python
if v(k) is not None:
    near(k, v(k), spec[gk])
```

A missing row is therefore a **skip**, not a failure. A referee deleted five rows from one table
and the verifier reported **106/106 passed**.

Rewriting 38 call sites would fix the instances. Asserting **how many checks ran** fixes the
category:

```python
floor = g["integrity"]["min_domain_checks"]
n_ran = len(checks)
record(n_ran >= floor, "number of domain checks that ran", n_ran, f">= {floor}")
```

Demonstrated: deleting two rows from `multidonor_specificity.csv` silently disables **7 checks**
in a section nobody had touched, and the floor is the only thing that notices. It has since
caught its own author three times.

## 3.2 A check that an attacker could switch off by deleting a file

`recompute.py` is the one component that *proves* rather than reproduces: it rebuilds published
AUROCs from committed per-example scores. Its strongest block opened with

```python
if reh.exists() and REHEARSAL.exists():
```

so **deleting the evidence directory made the check pass silently**. Only zeroing the scores
failed it. Now both absences raise.

## 3.3 A duplicated number nobody compared

`verify.py` gated R1 on `cost_of_matching.csv`. That table is a **join of two others**,
`rehearsal_binding_gc.csv` and `rehearsal_binding_dinuc.csv`, and `grep rehearsal_binding
scripts/verify.py` returned **zero hits**. Nothing asserted the copies agreed.

Permuting `rehearsal_binding_dinuc.csv` against its dataset labels passed **166/166** while
turning "larger in 88/94" into 67/94, "94/94 fall" into 80/94, and the Wilcoxon p from 3.8e-17
into 1.5e-12, the last of which violates the project's *own* threshold, in a check that never
ran on the corrupted table.

Eight cross-table assertions at exact equality now close it. The attack fails at max|diff| 0.282.

## 3.4 Evidence tables that were load-bearing and ungated

Five summaries are regenerated by `run.sh` with `--from-cache`, which reads a committed
`*_per_dataset.csv` instead of redoing hours of refits. `grep per_dataset scripts/verify.py`
returned nothing.

Zeroing every per-arm gain column in `k_sweep_per_dataset.csv` and every AUROC in
`strand_placebo_per_dataset.csv`, then rebuilding, reproduced both summaries **bit-for-bit** and
passed. So the claim "the headline was rebuilt from raw sequence to within 1.2e-06" reproduced
from a table in which all 188 per-arm gains were zero.

The fix asserts the arithmetic linking evidence to summary, and cross-checks the k=4 gains
against the published rehearsal values. The attack now fails on 6 assertions.

**The teaching point running through all four.** A gate is worth exactly what an attack against
it proves. Every one of these passed inspection and failed a five-line script. Before claiming a
check works, corrupt its input and watch it fail.

---

# Part 4. Auditing the manuscript, and a checker that checked nothing

## Why

R1's primary contrast, +0.0397, appeared in the manuscript, in **no committed table**, under
**no golden key**, produced by **no script that could be found**. It could not be recomputed and
it could not fail, and it survived six rounds of adversarial review in that state, because a
reviewer reads a number rather than goes looking for it.

`scripts/audit_manuscript.py` extracts every number quoted to three or more decimals and reports
the ones traceable to no committed table and no golden key.

## The first version was nearly worthless

It pooled every numeric cell of every table, including two with 66,010 rows. With that much
data, the four-decimal grid over [0.5, 1.0] came out **73.9% saturated**: essentially any
AUROC-like number matched something by coincidence. A fabricated `0.9427` injected into the
manuscript **passed cleanly** while the tool reported a reassuring 9 orphans.

It also wrote its report into `results/tables/` and then read it back on the next run, finding
every orphan inside its own output and printing **zero**.

## The fix

The haystack is now golden keys, cells of genuinely small summary tables, and column aggregates
,  the operations a manuscript actually performs. Saturation dropped **73.9% -> 6.3%**. Matching
happens at the token's own precision, so "1.036" matches a stored 1.0357967. The script prints
its own false-negative rate on every run, including the uncomfortable part: the three-decimal
grid is 44% occupied, so a number quoted to three decimals is close to a coin flip, and
**percentages written to one decimal are invisible to it entirely**.

    orphans: 45 -> 2

and the two survivors are the algebraic identity that retracted an earlier claim, i.e. evidence
*about* a retraction rather than a claim.

**The teaching point.** A checker that will not state its own false-negative rate is asking to
be over-trusted. Measure it, print it, and let the reader discount accordingly.

---

# Part 5. Choosing the paper

Four candidates went to a four-member council: **A**, an RBP benchmarking note on R1; **B**, a
methods paper on three ways a benchmark delta misleads; **C**, a negative-results calibration
paper; **D**, ship the repo and stop.

**B was tempting and I argued for it.** On this data, two of the three "traps" *invert the
conclusion*: the log-odds contrast flips sign, and a near-zero attenuation reads as independence
when it means the opposite. That felt like a paper.

**B is dead.** Each trap has a published owner, and two remedies are better than mine:

| trap | owner |
|---|---|
| AUROC scale compression | Pencina, D'Agostino, Pencina, Janssens & Greenland 2012, *Am J Epidemiol* 176(6):473-481 |
| coefficients not comparable across fits | Karlson, Holm & Breen 2012, *Sociological Methodology* 42:286-313, with a significance test mine lacks |
| attenuation against a calibrated null | Schuster, Twisk, ter Riet, Heymans & Rijnhart 2021, *BMC Med Res Methodol* 21:136 |

and Whalen, Schreiber, Noble & Pollard 2022, *Nat Rev Genet* 23:169-181 occupies the
"bundle of genomics ML pitfalls" slot at a venue nobody here can compete with.

**The correction to my own argument was the useful part.** Demonstrating that these traps
*invert conclusions* establishes **severity**, not **novelty**, and severity was never in
doubt. I was also quoting "2 in 3" as a base rate when n = 3, chosen post hoc from my own
retractions.

**The teaching point.** "I found three instances of X" is not a contribution if X is published.
Check who owns the remedy before building a paper on it.

---

# Part 6. The strand artifact: a control, a placebo, and a retraction inside the control

## The bug

`annotation.py:126` says so itself:

> Strand is deliberately dropped. A window's strand comes from its peak, so the region's own
> strand is never needed for classification or negative matching.

True for positives. False for negatives. `negatives.py:328` therefore writes

```python
"strand": p["strand"], "seq_rna": to_rna(dna, p["strand"])
```

giving each negative the **positive's** strand. Only **47.4%** of negatives demonstrably sit on
the strand their own gene is transcribed from, while every positive is true sense RNA. Direction
is a cue that separates the classes for a reason that is not binding.

## Why the cheap test could not settle it

The first attempt regressed each dataset's contrast on its sense fraction: rho = **-0.24
[-0.54, +0.11]**. Underpowered by construction, `frac_sense` spans only 0.433 to 0.615, so a
*between*-dataset regression is being used against a bias present *within* every dataset, and
the regressor is not exogenous (it correlates +0.427 with GC-arm AUROC, so it proxies region
mix). A linear extrapolation to fully-correct stranding gave **-0.0458 [-0.189, +0.072]**: an
interval containing both the published contrast and its negation.

Worth recording: a **$25 regeneration had been cancelled** on the evidence that the contrast
*grew* on sense-only negatives, +0.2643 -> +0.2787. That number is the composition-share
contrast, retracted since as an algebraic identity. **The cancellation rested on withdrawn
evidence and nobody noticed until it was checked.**

## The pre-registration

Criteria were written into `docs/61` and committed **before** the experiment ran: sign retained,
interval excluding zero, and at least 60% of the point estimate surviving. Restriction moves the
sense fraction to 1.0 by construction, which is full leverage.

## Why restriction alone lies, and the placebo

Keeping only pairs whose negative is unambiguously sense discards **57%** of pairs. A
256-feature k-mer model loses more from that than a 19-feature composition baseline does, **in
both arms**. So the contrast shrinks whether or not strand matters. Restriction alone reports
**-0.0091**, and a naive reading attributes all of it to strand.

Dropping the *same number* of pairs at random reproduces **-0.0032** of that shrinkage with no
strand involved.

## The placebo itself was not exchangeable

Uniform random dropping is the wrong counterfactual for a restriction that is not random.
Sense-kept negatives are **more intronic** (0.4335 against 0.4024 dropped) while GC is balanced
(0.5314 against 0.5332). Matching the placebo on the retained **region** marginals accounts for
a further **-0.0024** [-0.0048, -0.0003].

Region was not enough either. Retention *requires exactly one overlapping gene strand*, so it
selects against multi-gene loci by construction; gene density differed at a standardised mean
difference of **-0.303**. The final design stratifies on **region x gene-density quartile**.

## The retraction inside the control

With the unstratified placebo the excess was -0.0059 with an interval clear of zero, and the
section said **"the artifact is real"**. Matched on region the interval includes zero. That
sentence was withdrawn and only a **bound** is claimed.

**The teaching point, and it is the most transferable thing in this document.** A "restrict and
recompute" control needs a matched placebo, because restriction changes the sample size and the
sample composition at once. And the placebo needs stratifying on whatever the restriction
correlates with, which you must *measure*, not assume. Each layer here was found only by
someone checking the layer above.

## An argument that needs no placebo at all

`scripts/strand_asymmetry.py`, coordinates only, no model fitting. On the 40 datasets canonical
in both arms, GC-matched negatives are **42.8%** sense against **47.4%** in the dinucleotide
arm: paired difference **+0.047**, higher in **37/40**, Wilcoxon p = 5.0e-09.

The spurious cue is **stronger in the arm with the smaller gain**, so it works *against* the
reported direction. **The contrast is conservative, not inflated.** This is checkable in
seconds and depends on none of the placebo machinery.

---

# Part 7. What the magnitude is worth, and the paper's only biology

R1 concedes that the *sign* of the contrast is implied by the design: the composition baseline
spans the 15-df dinucleotide simplex, and that is exactly what the dinucleotide matcher
controls. Only the magnitude is informative. Three sections answer "what is the magnitude
worth".

## R1d, replication and efficiency

Fifteen proteins were assayed in **both** HepG2 and K562, separate eCLIP experiments with
separately drawn negatives, so one line is an out-of-sample prediction of the other. The
contrast replicates at **r = +0.909** [+0.812, +0.972].

The design guarantees the sign in each line independently and guarantees **nothing** about the
magnitudes agreeing, which is what makes this the answer to the design-implied-sign objection.

Two things were claimed too strongly and are now corrected:

- **"Replicates better than either arm"** was asserted *and gated*. A paired protein bootstrap
  gives r_contrast - max(r_arm) = **+0.113 [-0.082, +0.465]**, P(<=0) = 0.21 at n = 15. The gate
  is now inverted: the interval must **straddle** zero.
- **Most of r = 0.909 is protein signal strength.** The contrast correlates with total nested
  gain at **+0.952**. Partialling that out, replication falls to **+0.332** [-0.116, +0.690],
  p = 0.227. R1f already carried exactly this control and R1d did not, while being offered as
  the answer to the paper's biggest concession. Now gated as a limitation.

Efficiency: z = gain/SE rises **1.31x**, higher in 83/94. Because z grows as sqrt(n) that
converts into sample size, but carefully. The **ratio of means** gives 58% of the windows, the
**median dataset** 59%, and the **mean over datasets 96%**, with **14/94 needing more**. All
four forms are gated so the favourable one cannot be quoted alone.

## R1e, rebuilt from sequence, and independent of k

Every other R1 number is *read* from a table the analysis pass wrote. `scripts/k_sweep.py`
refits composition and k-mer models from `dataset.tsv` sequence for all 94 datasets in both
arms:

    rebuilt +0.0397   committed +0.0397   difference 1.2e-06   zero datasets skipped

| k | GC gain | dinuc gain | contrast | dinuc larger |
|---|---|---|---|---|
| 3 | 0.0168 | 0.0479 | +0.0311 | 90/94 |
| 4 | 0.0265 | 0.0662 | +0.0397 | 88/94 |
| 5 | 0.0277 | 0.0674 | +0.0397 | 89/94 |
| 6 | 0.0245 | 0.0600 | +0.0355 | 91/94 |

Positive at every k, and at every k in **82/94** datasets individually.

**This also settled an embarrassment quantitatively.** The manuscript called the model a
**5-mer** in four places. Both rehearsal tables record **k = 4** on all 189 rows; the "5" came
from `config/params.yaml` `cv: k: 5`, which is the **cross-validation fold count**. It survived
150 gated assertions because every one of them checks a *value* and none checked *what produced
it*. The table shows what the error would have cost: **k=5 minus k=4 = +0.0001, p = 0.84.**

## R1f, the only biology

| dominant region | n | contrast | composition alone |
|---|---|---|---|
| CDS | 24 | **+0.0635** | 0.5765 |
| 3'UTR | 17 | +0.0328 | 0.5971 |
| intron | 49 | **+0.0316** | 0.6656 |

CDS minus intron **+0.0319** [+0.0190, +0.0453], p = 1.5e-05.

The **mechanism was a prediction and it held**: intronic sites are compositionally distinctive
(polypyrimidine tracts, U-rich stretches), so composition alone should already discriminate them
and leave the protocol less to expose. Composition-only AUROC is **0.6656** intron against
**0.5765** CDS, p = 2.7e-08. That ordering is gated, so if it ever inverts the explanation fails
loudly.

**The limitation is gated beside it.** Partialling out total nested gain removes the association
(+0.082, p = 0.435). Region indexes *how much* non-compositional signal exists; it is not an
independent mechanism.

---

# Part 8. The cut

Six sections were removed. Their code and tables remain in the repository, gated and
reproducible, replaced in the manuscript by one table saying what was rejected and why.

The cut did more than tidy. **All 30 remaining manuscript orphans lay in cut sections; none was
ever in R1.** Removing them, and pruning the five limitations that belonged to those sections,
took orphans to **2**.

Three defects in R4d were fixed **even though R4d was being cut**, because the script ships:

1. The null drew a **normal** covariate when real phyloP is skewed (+1.04, sd 3.19, range -20 to
   +10). Non-collapsibility depends on the covariate's whole distribution, not its variance.
   Drawing through a Gaussian copula onto the empirical marginal moved the null from **-68% to
   -87%**.
2. The script **asserted** that a closed-form approximation "must agree" with the simulation and
   gated it at 0.15. That approximation is a small-sigma probit form; at c = 2.12 it is far
   outside its range and **anti-conservative**, so the gate was checking the correct simulation
   against a known-invalid reference and would have fired on the right answer.
3. Seeds 12 -> 40, against a seed-to-seed sd of 0.054.

**The teaching point.** Do not ship a script that asserts something you know is false, even if
no published claim depends on it. Someone will run it.

---

# Part 9. The state, and what is left

**241 gated assertions, 594 tests, 9 figures, 42 tables, 2 manuscript orphans, ~$38 spent.**

The paper is R1 plus five controls: R1b (scale), R1c (strand), R1d (replication and efficiency),
R1e (rebuild and k-robustness), R1f (region heterogeneity), with R2 retained as a methods table
citing Horlacher as replication.

Three writing tasks remain, all deferred deliberately: converting the audit document into a
manuscript with the usual sections, the README, and scrubbing a billing identifier before
publication. The branch `r1-scale-check-and-corrected-refit` is unmerged to `master`.

---

# Part 10. What to actually learn from this

These are ordered by how much they cost to learn.

1. **A control is not evidence until something has tried to break the control.** Every serious
   finding in this period was in work already reported as complete, and every one was found by
   someone running code, not by re-reading.

2. **Restriction is not random.** Any "restrict and recompute" test needs a matched placebo, and
   the placebo needs stratifying on whatever the restriction correlates with. Measure that; do
   not assume it.

3. **A duplicated number is worth exactly the assertion that the copies match.** Without it,
   permuting a source table passed 166/166.

4. **Never compare a subset rebuild against a full-panel published mean.** This panel-mixing
   error appeared three separate times in three different analyses.

5. **Gates written `if value is not None` disappear rather than fail.** Assert how many checks
   ran; it fixes the category rather than the instances.

6. **Check whether corroboration is independent before calling it corroboration.** Two
   derivations sharing inputs and algebra are one derivation.

7. **When an estimate needs a modelling choice, report the family.** A confidence interval
   conditional on a choice is not a confidence interval for the quantity.

8. **A checker should print its own false-negative rate.** Ours is 6.3% at four decimals and
   ~44% at three, and percentages are invisible to it. Saying so is what makes the 2 meaningful.

9. **Every assertion in this project checks a value; none checked what produced it.** That is
   how a 4-mer was called a 5-mer in four places and passed 150/150.

10. **Severity is not novelty.** Demonstrating that a known pitfall inverts a conclusion is
    useful, and is not a paper if the remedy is already published.
