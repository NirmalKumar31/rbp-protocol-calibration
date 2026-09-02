# 60. The paper: every claim, every number, every limitation

**Manuscript source of truth.** Written 2026-08-27 after seven rounds of adversarial review.
If a number in the manuscript disagrees with this file, this file is wrong and must be fixed,
not the manuscript quietly edited.

**What the harness does and does not establish.** `scripts/verify.py` runs **378 numeric
assertions** against `config/golden.yaml`, and the number of assertions that ran is itself
asserted, so a gate cannot silently skip. Most of those are **regression gates**: they detect
that a number changed, not that it was correctly derived. Two are stronger.
`scripts/recompute.py` rebuilds 285 published AUROCs from committed per-example scores at
max|difference| 2.2e-16 and fails if that evidence is zeroed or deleted, and `scripts/k_sweep.py`
rebuilds the headline contrast from raw sequence to within 1.2e-06. Call the rest a
reproducibility harness, not verification.

`scripts/audit_manuscript.py` reports every number in this document traceable to no committed
table and no golden key. It stands at **2**, both of them the algebraic identity that retracted
the composition-share claim, i.e. evidence about a retraction rather than a claim. **Its known
blind spot:** it only checks values written to three or more decimals, because the coarser grids
are too saturated to discriminate, so a percentage written as "55.2%" is invisible to it and has
to be checked by hand.

**Working title:** *What a sequence model adds over composition is set by the benchmark's
headroom: a three-protocol calibration across 94 ENCODE eCLIP datasets*

*Framing note.* An earlier title was "There is no protocol-independent measure of what a
sequence model contributes". Two independent reviewers recommended against it for the same
reason: a universal negative invites a reader to go find a counterexample, and three protocols
cannot support a universal quantifier. The calibration framing states the same evidence, keeps
the negative as the closing line rather than the title, and is the only framing in which the
mechanism (Spearman −0.60 against the baseline) is an asset rather than a confession.

**Venue:** bioRxiv immediately; then *NAR Genomics & Bioinformatics* (Methods/benchmarking).

---

## The one-sentence claim

> **The composition baseline your negative set leaves behind is most of what determines the
> measurable contribution of a sequence model.** Across 94 ENCODE eCLIP datasets, holding the model, the positives,
> the folds and the estimator fixed and changing only how the negatives were built, the nested
> contribution of a 4-mer over a 19-feature composition baseline measures **+0.0663**,
> **+0.0265** or **+0.0122 AUROC** across three protocols -- a **5.4-fold range** [4.4, 6.6] --
> and falls monotonically as the baseline rises. **Given the baseline, knowing which protocol
> produced it adds 1.0% of variance; given the protocol, knowing the baseline adds 11.0%.**
> Where the two disagree, the baseline wins: in the 27 of 94 datasets where the bias-aware
> protocol *lowers* the baseline relative to GC matching, its contribution deficit reverses
> (−0.0212 → **+0.0028**), and matched on baseline the headline dinucleotide-versus-GC contrast
> falls from +0.0398 to **−0.0087 [−0.0265, +0.0122]**. **No rescaling recovers a
> *transportable* protocol-free quantity.** Over eight standard monotone transforms the range
> floors at **2.00x** [1.67, 2.46], far outside the equal-means null (95th percentile 1.18);
> an exponent that equalises our three protocols does exist (p = 1.544, range 1.00x) but it is
> a property of the benchmark, not of the quantity -- **the exponent fitted on an independent
> benchmark is 3.65, and ours leaves that benchmark at 2.34x against the 2.38x it started
> from.** The practical consequence is to report the composition-only AUROC under the same
> protocol alongside every headline AUROC, and never to compare contributions across
> protocols.

**The one qualification the claim must carry, and it is not optional.** "Most" is not "all".
The 1.0% is an average over protocol pairs that behave differently: **0.05%** for GC versus
dinucleotide, **3.40%** for GC versus neg2, and at matched baseline neg2 still costs
**−0.0081 [−0.0130, −0.0036]**. The baseline is essentially the whole story *within* the
composition-matched family and not across families, because the baseline gradient itself is a
property of composition-matched negatives (−0.545, −0.462) and vanishes for other-RBPs'-sites
negatives (−0.122, n.s.). An earlier version of this claim said the protocol label carries
essentially no information; that is **withdrawn**, and R1n carries the corrected form.

**Deliberately NOT in the one-sentence claim, after two rounds of statistical review.** The
decomposition into "compression" and "protocol effect" is a supporting sensitivity analysis,
not a headline quantity: R1h shows it is not identified for any model. The raw contrast and the
multiplier need no transplant, no link and no transportability assumption, and they are what
survived every attack.

---

## R1: the protocol effect (primary result)

`cost_of_matching.csv`, n = 94 paired datasets. Figures **f0** (the panel) and **f1**.

| quantity | GC-matched | dinuc-matched |
|---|---|---|
| composition alone (19 features) | 0.7827 | 0.6280 |
| 4-mer score alone | 0.7981 | 0.6886 |
| composition + 4-mer score | 0.8092 | 0.6941 |
| **nested contribution of the score** | **+0.0265** [+0.0210, +0.0324] | **+0.0662** [+0.0559, +0.0765] |
| datasets where the score adds significantly | **80/94** (72/94 at the measured design effect) | 82/95 (78/95) |
| datasets where composition ≥ the k-mer model | **29/94** | n/a |

An earlier draft of this table put 0.7981 and 0.6886 on the "composition + score" row. Those
are the *standalone* 4-mer AUROCs. The nested values are 0.8092 and 0.6941, and the reported
contributions were always the nested ones, so the claim did not change; the row label was
wrong and is corrected here.

- change in apparent AUROC under dinucleotide matching: **−0.1095**, **94/94 datasets fall**, paired Wilcoxon p < 1e-15
- contrast in nested contribution: **+0.0397** [+0.0336, +0.0458]
- **effect modification by dataset size is real: rho = +0.307, p = 0.0026** on the reported
  nested contrast (Figure **f5**). The paired design still rules out size *confounding* (both arms use the
  same datasets), but the effect is larger in larger datasets and that must be printed. (An
  earlier draft printed rho = 0.141, p = 0.175, which is this correlation computed on the
  RETRACTED composition-share difference, not on the quantity the paper claims.)

**The circularity caveat, which must appear in the same paragraph.** The composition baseline's
19 features are mononucleotide + dinucleotide counts + entropy. **GC is not one of them**; it
is the spanned combination C+G. Mononucleotide counts are marginals of dinucleotide counts and
contribute no independent degrees of freedom, so the baseline's frequency space is the 16-cell
dinucleotide simplex: **15 df**. GC is one linear functional of it. So the GC arm controls **1
of 15** and the dinucleotide arm controls **15 of 15**. (An earlier draft said 1-of-18 and
18-of-18 while simultaneously conceding mono is spanned by dinuc; both cannot hold.) The *sign* of the +0.0397 contrast is therefore
implied by the design; only its *magnitude* is informative. Say exactly this. Do not lead with
the 94.8% / 67.8% share framing; see the retraction below.

## R1b: the contrast is not an artefact of the AUROC scale

`scale_check.csv`, `scripts/scale_check.py`, 13 gated checks. Figure **f8**.

**The objection.** R1 compares two nested AUROC gains whose baselines differ a lot: composition
alone reaches 0.7827 in the GC arm and 0.6280 in the dinucleotide arm. AUROC is bounded at 1,
so a fixed increment of real discriminative signal buys a smaller AUROC increment when the
baseline is already high. The measured compression factor between these two baselines is
**1.51x**. The objection is therefore not vague: scale compression predicts the observed
direction with no protocol effect at all, and it has to be answered with a number.

Somers' D does not answer it. D = 2·AUROC − 1, so the nested gain on the D scale is exactly
twice the same quantity and the contrast merely doubles. It is a linear rescaling of the thing
in question.

**The decomposition.** Under a binormal model d′ = √2·Φ⁻¹(AUROC) is linear in signal and
unbounded. Transplanting the GC arm's own d′ increment onto the dinucleotide arm's baseline
predicts what that arm would show if the protocol moved the baseline and changed nothing else.

| transplant direction and link | compression | **protocol effect** |
|---|---|---|
| GC increment onto the dinuc baseline, probit | +0.0083 [+0.0061, +0.0109] | **+0.0313** [+0.0267, +0.0363] |
| dinuc increment onto the GC baseline, probit | +0.0182 [+0.0142, +0.0223] | **+0.0215** [+0.0180, +0.0252] |
| GC onto dinuc, logit link | | +0.0288 [+0.0243, +0.0335] |
| dinuc onto GC, logit link | | **+0.0188** [+0.0156, +0.0223] |

**Report the range, not the point.** The same logic licenses transporting either arm's
increment, and under either link, and the four disagree: compression accounts for 21% of the
contrast in the most favourable member and 46% in the least. That choice moves the estimate
further than any single interval is wide, so quoting **+0.0313** alone would be question-begging.
Each row's interval is a legitimate bootstrap interval *conditional on that transport choice*;
none of them has coverage for "the protocol effect" as such.

> **SUPERSEDED BY R1h, AND THE CORRECTION IS LARGE.** This section previously concluded "the
> protocol effect is therefore **+0.0188 to +0.0313**, and every member of the family keeps the
> sign". Both halves were too strong. Four members is an arbitrary truncation of the link
> family: over six ROC-motivated links x two directions the honest range is **+0.0127 to
> +0.0506**, and the odds scale -- a member of the same power family, excluded only because it
> has no ROC derivation -- gives **−0.0036**. More seriously, the transplant's identifying
> assumption is false: the d′ increment is not baseline-invariant, and once that is corrected
> **no model's protocol effect is identified**. Read this table as a description of what the
> four originally-chosen members give, not as a bound. R1h is the current claim.

The protocol effect is positive in **87/94** datasets and its interval excludes zero. The same
comparison computed directly on the unbounded d′ scale gives **+0.1290** [+0.1091, +0.1499],
larger in 87/94, which would be zero if compression were the whole story. Corrected for
compression, the GC-matched arm reports **47.4%** less of the model's measurable contribution
than the dinucleotide-matched arm does. The uncorrected figure is 60%, and an earlier draft
said "two-thirds", which was wrong on both counts. Note the phrasing: neither arm is the truth,
so the figure describes a difference between two protocols and not something one of them hides.

**The reversal, which must be reported.** The same comparison on a third scale, the Firth
coefficient of the standardised score in the nested fit, **changes sign**: +1.063 in the GC arm
against +0.686 in the dinucleotide arm, dinucleotide larger in only 11/94. A result whose sign
depends on the scale is not a result unless the reversal itself has a diagnosis, and it does. A
logistic coefficient is identified only against the latent residual scale, so coefficients from
two fits with different total signal are not comparable (Mood 2010, *Eur Sociol Rev*; the same
non-collapsibility that voided this paper's earlier conservation analysis, R4c). The GC task
carries **1.79x** the total signal of the dinucleotide task. The fingerprint is direct: across
the 94 datasets the between-arm coefficient gap tracks each task's **total** signal at
Spearman **+0.520** (p = 8e-08) and tracks the **incremental value** it is supposed to measure
at **+0.065** (p = 0.53). Dividing each coefficient by its own fit's total signal restores the
sign, **+0.1154** [+0.0715, +0.1597], larger in 68/94.

*A correction to an earlier draft of this section.* It also reported a "scale-null residual" of
+0.1156 and presented it as a second, independent answer. It is not independent: the residual is
algebraically identical to the normalised contrast above under a positive weight (verified at
max|difference| 5.0e-16, sign agreement 94/94). And the null it is measured against, that the
coefficient scales as the first power of total signal, is a choice rather than a derivation:
under a square-root scaling the same residual is **−0.1721** and under a 3/2 power it is
**+0.5278**, so its sign is a free exponent. Both numbers remain in the committed tables and
nothing is claimed from them. The load-bearing answer to the reversal is the fingerprint, which
does not depend on any assumed scaling.

So two of three scales agree with R1 and the third is measuring task difficulty rather than
incremental value. The reversal and its fingerprint are gated, so if the fingerprint ever
inverts the verifier fails and the reversal becomes real evidence against R1.

**Provenance note, stated because it matters.** Before this analysis the +0.0397 contrast
appeared in the manuscript, in no committed table, and under no golden key. It could not be
recomputed and it could not fail. It is now derived by `scripts/scale_check.py` from the two
rehearsal arms and asserted by `r1_scale_check`. Its interval, recomputed with a paired
dataset-level bootstrap (2000 draws, seed 0), is [+0.0336, +0.0458]; the manuscript previously
printed [+0.0334, +0.0461] from an unrecorded computation.

## R1c: the strand artifact, measured against a matched placebo

`strand_placebo.csv`, `scripts/strand_placebo.py`, 16 gated checks, n = 40 datasets.
Figure **f3**.

**The wound.** `annotation.py:126` drops region strand by design ("a window's strand comes from
its peak"), which is right for positives and wrong for negatives, so `negatives.py:328` gives
each negative the *positive's* strand. Only **55.2%** of negatives sit on the strand their own
gene is transcribed from, counting unambiguous windows only (**47.4%** of all negatives); positives are 100% true sense. Direction is therefore a cue that
separates the classes for a reason that is not binding, and it inflates absolute AUROCs in both
arms.

**Why the obvious control lies, twice.** Keeping only pairs whose negative is unambiguously
sense discards 57% of pairs. A 256-feature 4-mer model loses more from that than a 19-feature
composition baseline does, in both arms, so the contrast shrinks whether or not strand matters.
And the pairs that survive are not a random half: their negatives are more intronic (**0.4335**
against **0.4024** dropped in the GC arm, 0.4413 against 0.4003 in the dinucleotide arm) while
GC is balanced (**0.5314** against **0.5332**). So the restriction reweights the task by locus
type, and a uniform-random placebo is the wrong counterfactual for it. That is why the placebo
is stratified on region and not on GC (`strand_asymmetry.csv`).

| | contrast, n = 40 |
|---|---|
| full data | **+0.0378** [+0.0288, +0.0478] |
| sense-only pairs (43% retained) | +0.0287 [+0.0201, +0.0383] |
| placebo, same n dropped at random | +0.0338 [+0.0250, +0.0436] |
| placebo, matched on region x gene-density | +0.0342 [+0.0256, +0.0436] |

Decomposing the −0.0091 that restriction alone reports: **−0.0040** [−0.0070, −0.0019] is the
cost of discarding pairs at all, and the remainder, **−0.0055** [−0.0089, −0.0022], is strand.
The strand-corrected contrast is **+0.0322** [+0.0237, +0.0421] and **85.4%** of the effect
survives. The artifact is **real** and small: its interval excludes zero.

**Two corrections, and the order they happened in matters.** An earlier version of this section
matched the placebo on region alone with five seeds per dataset per arm. It reported an excess
of −0.0036 with an interval of [−0.0071, +0.0001] that touched zero, and withdrew the claim that
the artifact is real. Two things were wrong with that. First, region alone was not the right
stratum: retention *requires exactly one overlapping gene strand*, so it selects against
multi-gene loci by construction: a sense-kept negative overlaps **1.10** annotated genes on
average against **1.37** for a dropped one. Second, and larger, five placebo seeds
left roughly a sixth of the between-dataset variance as Monte Carlo noise, which inflated the
interval. At twenty seeds and region-by-gene-density strata the excess is −0.0055 with an
interval clear of zero, and the claim is restored.

A third thing follows and is worth stating because it cuts against the earlier reasoning: with
adequate seeds, **the stratification barely matters**. The stratified and unstratified placebos
now differ by **+0.0004** [−0.0010, +0.0019]. The −0.0024 "locus mix" reported earlier was
largely seed noise, not locus mix. The placebo design was less load-bearing than one noisy run
made it look, and the paper gates that gap as a bound rather than as a finding.

**The artifact cannot manufacture the contrast, only shrink it, and this does not depend on the
placebo.** On the 40 datasets whose window tables are canonical in both arms, GC-matched
negatives are **42.8%** sense (sd 0.018) against **47.4%** in the dinucleotide arm: a paired
difference of **+0.047**, higher in **37/40**, Wilcoxon p = 5.0e-09 (`strand_asymmetry.csv`;
percentages count ambiguous-strand windows in the denominator, the same convention as the
47.4% quoted above). The spurious directional cue is therefore *stronger in the arm with the
smaller nested gain*. It works against the reported direction, so **+0.0397 is conservative**
rather than inflated. This argument uses coordinates only, no model fitting, and is independent
of the placebo experiment.

**Pre-registered.** The criteria were committed before the experiment ran: sign retained,
interval excluding zero, at least 60% of the point estimate surviving. All three hold on the
primary (stratified) estimator. Every dataset used is gated on reproducing its own published
row first; the dinucleotide arm's canonical window tables were fetched for this purpose because
the local copy was a different draw and reproduced only 13 of 40.

**A weaker version is retained in `strand_contrast.py` and should not be relied on.** It
regresses the per-dataset contrast on each dataset's sense fraction. On the 40 audited datasets
it gives rho = −0.24 [−0.54, +0.11]; extended to all 94 it gives rho = −0.074, p = 0.478. It
has almost no power because `frac_sense` spans only 0.433 to 0.615, so a *between*-dataset
regression is being used against a bias present *within* every dataset, and `frac_sense` is not
exogenous. Restriction moves the sense fraction to 1.0 by construction. Both designs are kept
because the difference between them is the methodological point.

## R1d: what the magnitude is worth

`r1_robustness.csv`, `scripts/r1_robustness.py`, 11 gated checks. Figure **f4**.

R1 concedes that the *sign* of the contrast is implied by the design, so only the magnitude is
informative. This section is what the magnitude is worth. Neither result below is implied by
the design.

**It replicates out of sample.** Fifteen proteins were assayed in both HepG2 and K562. Those
are separate eCLIP experiments with separately drawn negatives, so the contrast measured in one
line is an out-of-sample prediction of the other. Pearson **r = +0.909** [+0.812, +0.972],
p = 2.6e-06; Spearman +0.932. Mean absolute difference between lines **0.0151**, against a
between-dataset standard deviation of 0.0318.

The contrast replicates **at least as well as either arm it is built from** (GC-arm gain alone
r = +0.518, dinucleotide-arm gain alone r = +0.813). Stated carefully, because the ordering
itself is not established: a paired protein bootstrap gives r_contrast − max(r_arm) =
**+0.113** [−0.082, +0.465], P(≤0) = 0.21 on n = 15. The claim is that differencing the two arms
does not destroy the signal, not that it improves on them.
The design guarantees the sign in each cell line independently and guarantees nothing about
whether the magnitudes agree, so this is information about magnitude alone.

**The limitation, and it is gated.** A protein's contrast is strongly related to how much
total non-compositional signal that protein has: the two correlate at **+0.952** (p = 5.0e-08).
So two cell lines agreeing on the contrast may be two cell lines agreeing on the protein's
signal strength rather than on a protocol-specific quantity. Partialling out the per-protein
mean total gain, the replication falls from +0.909 to **+0.332** [−0.116, +0.690], p = 0.227.
What is established at n = 15 is that the magnitude is a reproducible property of the protein;
what is *not* established is that it is reproducible **beyond** the protein's overall signal
strength. This is the same control R1f carries, and R1d needed it more, having been offered as
the answer to the design-implied-sign objection.

**It buys statistical efficiency, which is what a benchmark builder acts on.** The nested gain
rises 2.50x under dinucleotide matching, but its standard error rises only 2.18x, so the
per-dataset signal-to-noise z = gain / SE rises **1.31x** (mean 9.97 to 13.11, median 7.01 to
10.13), higher in **83/94**, paired Wilcoxon p = 3.1e-14. Because z grows as the square root of
sample size, that converts into sample size, but it must be reported carefully: the ratio of
means gives **58%** of the labelled windows, the **median dataset** needs **59%**, and the
**mean over datasets is 96%** because the advantage is concentrated rather than universal.
**14/94 datasets need *more* windows** under dinucleotide matching. The honest summary is that
the median dataset reaches the same confidence on roughly 60% of the data and a sixth of
datasets do worse. And composition alone beats the sequence model in **29/94** datasets under GC matching
against **14/94** under dinucleotide matching: the harder protocol makes the model look better
on that comparison, not worse.

**This is the practical recommendation, and it does not depend on either arm being "the
truth".** A benchmark built on dinucleotide-matched negatives reports a lower headline AUROC
and a more precisely measured, more reproducible, and roughly twice as large estimate of what
the model actually contributes.

## R1e: the contrast rebuilt from sequence, and independent of k

`k_sweep.csv`, `scripts/k_sweep.py`, 11 gated checks, all 94 datasets, both arms.
Figure **f6**.

**Rebuilt, not re-read.** Every other number in R1 is read from `rehearsal_binding_*.csv`,
which the analysis pass wrote, so its gate detects drift rather than error. This section refits
the composition baseline and the k-mer model from `dataset.tsv` sequence for all 94 datasets in
both arms and recomputes the contrast from scratch. Every dataset first passes a per-arm
reproduction gate at k = 4; **none was skipped**. The rebuilt contrast is **+0.0397** against a
committed **+0.0397**, absolute difference **1.2e-06**. With `recompute.py`, which rebuilds 285
published AUROCs from per-example scores, this is the second result in the study that is proved
rather than reproduced.

**Independent of the k-mer size.**

| k | GC-arm gain | dinuc-arm gain | contrast | dinuc larger |
|---|---|---|---|---|
| 3 | 0.0168 | 0.0479 | **+0.0311** [+0.0267, +0.0355] | 90/94 |
| 4 | 0.0265 | 0.0662 | **+0.0397** [+0.0336, +0.0458] | 88/94 |
| 5 | 0.0277 | 0.0674 | **+0.0397** [+0.0334, +0.0461] | 89/94 |
| 6 | 0.0245 | 0.0600 | **+0.0355** [+0.0291, +0.0422] | 91/94 |

The contrast is positive at every k, and positive at *every* k in **82/94** datasets
individually. So it is not an artifact of one arbitrary modelling choice.

**A disclosure this table also settles.** Four places in an earlier draft described the model as
a 5-mer. It is a 4-mer: both rehearsal tables record k = 4 on all 189 rows, and the "5" came
from `config/params.yaml` `cv: k: 5`, which is the cross-validation fold count. The error is
corrected throughout, and the table above shows what it would have cost: the k = 5 and k = 4
contrasts differ by **+0.0001**, paired Wilcoxon **p = 0.84**. Not one conclusion in this paper
turns on it. The reason it survived 150 gated numeric assertions is worth stating plainly,
because it generalises: every one of those assertions checked a *value*, and none checked
*which model produced it*. `k` is now gated.

## R1f: the effect is twice as large for coding-region binders

`region_heterogeneity.csv`, `scripts/region_heterogeneity.py`, 10 gated checks, n = 94.
Figure **f7**.

This is the only biological statement the study supports, and its limitation is printed beside
it rather than below it. Grouping datasets by the region their positive windows mostly occupy:

| dominant region | n | contrast | composition alone (dinuc arm) |
|---|---|---|---|
| CDS | 24 | **+0.0635** | 0.5765 |
| 3'UTR | 17 | +0.0328 | 0.5971 |
| intron | 49 | **+0.0316** | 0.6656 |

CDS minus intron is **+0.0319** [+0.0190, +0.0453], Mann-Whitney p = 1.5e-05, Kruskal-Wallis
across the three groups p = 2.2e-05. The difference between binding classes is as large as the
headline contrast itself.

**The mechanism is checkable and it checks out.** Intronic binding sites are compositionally
distinctive: polypyrimidine tracts, U-rich stretches, the low-complexity sequence around splice
signals. Composition alone should therefore already discriminate them, leaving the protocol
less to expose. It does: composition-only AUROC under dinucleotide matching is **0.6656** for
intron-dominant datasets against **0.5765** for CDS-dominant ones, p = 2.7e-08. Where
recognition is compositional the protocol change matters less; where it is higher-order, it
matters twice as much.

**The limitation, and it is gated so it cannot be dropped.** Region does not act independently
of effect size. The intronic fraction correlates with the contrast at Spearman −0.340
(p = 8.0e-04), but once the total nested gain is partialled out the association disappears
(+0.082, p = 0.435). So the honest claim is that **region indexes how much non-compositional
signal there is to expose**, not that region is a separate mechanism. Stated the second way it
is supported; stated the first way it would be wrong.

## R1g: the contrast is not a property of the model class, and survives at scale

*(This heading previously read "and it GROWS with capacity". **Withdrawn**: the ordering
reverses on the ratio scale, where the multiplier runs 3.08 / 3.51 / 2.38 and SpliceBERT is the
LOWEST, and R1h's specification grid has the CNN surviving most, not the largest model. What
survives is that all three model classes show the contrast, on all 94 datasets. See RETRACTED.)*

`deep_contrast.csv`, `deep_contrast_per_dataset.csv`, n = 94 paired datasets. Figure **f9**.

The GC arm was swept for two further architectures, 470 runs each: the 7,089-parameter CNN and
the 19.7M-parameter fully fine-tuned SpliceBERT. The dinucleotide arm's per-window scores were
already committed. Both arms use the same chromosome folds, the same seed, the same
hyperparameters, the same code path and full datasets; the only difference is how the negative
windows were chosen. The 4-mer is **refitted alongside them on the same rows** rather than read
from `rehearsal_binding_*`, so all three rungs are measured identically.

| model | size | nested, GC | nested, dinuc | contrast | positive | protocol effect |
|---|---|---|---|---|---|---|
| k-mer LR | 256 features | +0.0265 | +0.0663 | **+0.0398** [+0.0337, +0.0459] | 88/94 | +0.0188 to +0.0314 |
| CNN | 7,089 param | +0.0330 | +0.0860 | **+0.0530** [+0.0446, +0.0614] | 89/94 | +0.0213 to +0.0412 |
| SpliceBERT | 19.7M param | +0.0890 | +0.1754 | **+0.0864** [+0.0788, +0.0943] | 94/94 | +0.0253 to +0.0543 |

**The ladder is the claim, not mere survival.** The contrast grows monotonically with model
capacity, and because the three rungs are measured on the same 94 datasets the right statistic
is the paired difference:

| step | difference | larger in |
|---|---|---|
| CNN − k-mer | **+0.0132** [+0.0068, +0.0197] | 58/94 |
| SpliceBERT − CNN | **+0.0334** [+0.0289, +0.0380] | 90/94 |
| SpliceBERT − k-mer | **+0.0466** [+0.0416, +0.0517] | 94/94 |

Every step's interval excludes zero. **The CNN step is the weak rung** and is reported as such:
its mean difference is clear but it holds on only 58 of 94 datasets, against 94 of 94 for
SpliceBERT over the k-mer. Marginal intervals for the k-mer and the CNN overlap, which is
exactly why the paired form is the one reported.

**The compression correction is applied identically**, and it does NOT support a claim about
capacity. Each model's protocol effect is the same transplant family used in R1b, and every
member is positive for every model. The correction bites harder for the stronger model, as it
must: SpliceBERT's compression-only component is +0.0320 of its +0.0864, against +0.0083 of the
k-mer's +0.0398.

**An earlier version of this section claimed "the protocol effect still rises with capacity".
That claim is WITHDRAWN.** It was defended by putting two overlapping ranges side by side,
which is not a test; and when tested properly it rests entirely on an assumption R1h shows to
be false. Under the transplant corrected for baseline invariance, no model's protocol effect is
identified, and the model that survives the most specifications is the **CNN**, the middle
rung. Nothing about protocol sensitivity is monotone in parameter count on any scale. See R1h
for the specification grid.

### The ladder reverses on the ratio scale, and that is the diagnosis

R1b's rule is that **a result whose sign depends on the scale is not a result unless the
reversal itself has a diagnosis**. That rule was applied ruthlessly to the log-odds reversal and
then not applied to this section's own headline until a referee ran the ratio scale. Restricted
to the 77 datasets where every model has a positive gain in both arms:

| model | additive contrast | **protocol multiplier** |
|---|---|---|
| k-mer LR | +0.0457 | **3.08x** [2.62, 3.62] |
| CNN | +0.0619 | **3.51x** [2.87, 4.43] |
| SpliceBERT | +0.0935 | **2.38x** [2.12, 2.71] |

| paired step, log-ratio | difference | larger in |
|---|---|---|
| CNN − k-mer | +0.130 [−0.047, +0.337] | 42/77, n.s. |
| SpliceBERT − CNN | **−0.388** [−0.564, −0.246] | 14/77 |
| SpliceBERT − k-mer | **−0.258** [−0.369, −0.139] | 14/77 |

**The multiplier does not grow. It is smallest for the largest model**, and significantly so.
The additive ladder is real and the multiplicative one runs the other way, so "the contrast
grows with capacity" is a statement about the AUROC difference scale and must be said that way.

**The diagnosis.** The protocol multiplier is roughly invariant at **2.4x to 3.5x** across all
three model classes. Larger models have more total signal, so the same near-constant multiplier
produces a larger absolute difference; and the AUROC ceiling bites harder at SpliceBERT's higher
baseline, which is why its multiplier is the smallest of the three. Both facts are consequences
of one quantity, not two competing results.

**What this does to the paper.** The claim to make is the **invariance**, not the ladder:

> The negative-set protocol multiplies a model's measured contribution over composition by
> roughly 2.4x to 3.5x, and that multiplier holds across 79 RBPs, two cell lines and three
> model classes spanning 256 features to 19.7M parameters.

That answers "would this vanish for a real model?" decisively -- it is still 2.4x for a
fine-tuned transformer -- without the scale-dependent overclaim. The published 4-mer figure is
representative rather than special. It also unifies R1d and R1f, whose partial correlations both
die once total signal is removed: those are the same invariant multiplier applied to varying
amounts of signal, which is why removing the signal removes the association.

**Cross-check.** The refitted 4-mer reproduces R1's published contrast to **7.05e-05**, which is
what licenses comparing it against the deep models at all.

**Provenance, stated accurately after a referee checked it.** An earlier version of this
paragraph was wrong in two ways and both corrections weaken the CNN rung specifically.

*Where each arm ran.* SpliceBERT: dinucleotide on mixed A10G, GC on NVIDIA A10. **The CNN's
dinucleotide arm ran on x86 GCP Batch CPU, not a GPU** (`cloud/jobs/rendered/`,
`--device cpu`, `e2-standard-4`), while its GC arm ran on Apple MPS. So the CNN rung is a
CPU-versus-Metal comparison, which is the weakest provenance in the paper and sits on the
weakest rung. A referee's same-backend replication on 10 datasets moved that rung's contrast
UPWARD by roughly a quarter, so this does not overturn the ordering -- but a shift of that
size is comparable to the entire CNN − k-mer step and cannot be reported from a cross-backend
measurement without saying so. (Their figures are an external measurement, not reproduced in
this repository and therefore not gated here; the re-run below is what will settle it.)

*"Identical seed" was false.* `registry.build()` ran before `trainer.set_seed()`, so weights
were drawn from an unseeded RNG for the whole study and `torch.manual_seed(7)` governed only
dropout and batch order. All 945 fold-runs used an uncontrolled, unrecorded initialisation.
Fixed, and `tests/unit/test_seeded_init.py` now pins the ordering as well as the behaviour.
Measured cost: per-dataset SD of the nested contribution across independent training runs is
0.006 (CNN) to 0.010 (SpliceBERT), which induces ~0.001 on a panel mean over 94 datasets
against a reported CI half-width of ~0.008. **The panel means are unaffected; exact
reproducibility and the per-dataset counts are not.** Everything else claimed identical across
arms was checked and holds: seed value, epochs, batch size, learning rates, weight decay, fp32,
fold maps, and the training source files byte-for-byte.

*One further correction.* The dinucleotide arm has **no committed per-run metadata** -- only
scores -- so its epoch counts cannot be verified from what survives. And the sentence "the only
difference is how the negative windows were chosen" is false: the positive sets differ too
(pairwise Jaccard above 0.91 but below 1 on most datasets), because both matchers drop
positives they cannot match. Restricting
both arms to shared positives moves the contrast +0.0398 -> +0.0401, so it is immaterial, but
the sentence should say "almost the only difference".

The dinucleotide arm's windows drifted by **0.06%** (172 rows in 307,430) after its sweep was
scored, because negative matching is a stochastic search that was re-run; every model is
intersected to a common row set before anything is fitted, and the minimum coverage is gated.

## R1h: is the protocol effect identified? Two referee attacks, one of which lands

`protocol_identification.csv`, `protocol_baseline_slopes.csv`, n = 94.

R1b's transplant is the paper's answer to the compression objection, so it deserves the same
hostility everything else got. A statistical referee attacked it twice.

**Attack 1: "two directions, two links" is an arbitrary truncation, so the stated range is not
a range.** Correct, and the paper's numbers were too narrow. Under six standard monotone links
(probit, logit, arcsine, complementary log-log, −log(1−a), log a), both directions, twelve
members:

| model | published range | **honest range over 12 members** | sign holds |
|---|---|---|---|
| k-mer | +0.0188 to +0.0314 | **+0.0127 to +0.0506** | yes |
| CNN | +0.0213 to +0.0412 | **+0.0148 to +0.0589** | yes |
| SpliceBERT | +0.0253 to +0.0543 | **+0.0159 to +0.1022** | yes |

The sign survives every member for every model, so R1b's qualitative claim stands. **The
quoted range must widen**, and this is what the paper now reports.

**Attack 2: the identifying assumption is false.** This one lands. Transplanting an increment
across baselines requires the increment to be baseline-invariant. Regressing the within-arm d′
increment on the within-arm d′ baseline:

| model | slope, GC arm | p | net of mechanism |
|---|---|---|---|
| k-mer | −0.0646 | 0.020 | −0.0643 |
| CNN | −0.0576 | 0.143 | −0.0572 |
| SpliceBERT | **−0.3417** | **3.4e-07** | **−0.3407** |

Part of any such slope is mechanical: the baseline appears on both sides, so estimation noise
forces a negative slope. Computing that term correctly needs the **covariance** between the two
AUROC estimates, which is 0.38 to 0.86 here because both models are fitted on the same rows.
A first version of this analysis set it to zero, over-subtracted by 5x to 30x, and biased the
protocol effect high. The tell was that the mechanical slope came out identical across all
three models within an arm; it is now gated against exactly that.

**The transplant has three defensible slope estimators and they disagree.** The slope can be
taken from the source arm, the target arm, or pooled, and the two arms differ by 2x–3x
(k-mer −0.064 vs −0.207; SpliceBERT −0.341 vs −0.454; the CNN changes sign). A correctly
specified linear-in-baseline model would give the same slope in both arms. It does not, so the
adjustment is itself misspecified and what follows is a **sensitivity band, not a correction**:

| model | span across 6 specifications | keeps its sign in |
|---|---|---|
| k-mer | −0.0044 to +0.0200 | **4/6** |
| **CNN** | **+0.0190 to +0.0588** | **6/6** |
| SpliceBERT | −0.0103 to +0.0044 | **1/6** |

**What this costs the paper, stated at full strength.** **No model's protocol effect is
identified.** An earlier version of this section said "the k-mer's protocol effect survives
it"; that is true only under the source-arm slope, which is the most favourable of three
defensible choices, and it is withdrawn. Under the target arm's slope the k-mer gives
−0.0044 [−0.0104, +0.0015], the same failure conceded for SpliceBERT.

**And the ordering is refuted, not merely unsupported.** The model that survives every
specification is the **CNN** — the middle rung. Nothing about protocol sensitivity is monotone
in parameter count, on any scale: not the raw multiplier (R1g), not the apparent-AUROC drop,
and not the identified protocol effect.

**What it does not cost.** The raw contrast requires no transplant, no link and no
transportability assumption: it is a difference of two AUROC differences measured on the same
rows. **+0.0398, +0.0530 and +0.0864 stand**, as does the 2.4x–3.5x multiplier. Every attack in
two rounds of statistical review left them untouched. The identification problem is confined to
the decomposition, which is a supporting control and should be reported as a sensitivity band.

## R1i: every interval, recomputed with the clustering the panel actually has

`cluster_intervals.csv`, n = 94 datasets over **79 proteins**.

Fifteen proteins are assayed in both K562 and HepG2 and contribute two rows each. Every
bootstrap in this project resamples *datasets*, which treats those two rows as independent
evidence. They are not: the within-protein correlation of the contrast is **+0.924** for the
k-mer, +0.721 for the CNN, +0.761 for SpliceBERT. **Every interval printed above is therefore
narrower than the data support**, and the correct unit is the protein.

| quantity | by dataset (as published) | **by protein (correct)** | width |
|---|---|---|---|
| contrast, k-mer | [+0.0338, +0.0461] | [+0.0325, +0.0477] | 1.23x |
| contrast, CNN | [+0.0447, +0.0614] | [+0.0440, +0.0627] | 1.12x |
| contrast, SpliceBERT | [+0.0788, +0.0942] | [+0.0778, +0.0951] | 1.13x |
| step CNN − k-mer | [+0.0068, +0.0197] | [+0.0061, +0.0205] | 1.12x |
| step SpliceBERT − CNN | [+0.0289, +0.0380] | [+0.0287, +0.0383] | 1.05x |
| step SpliceBERT − k-mer | [+0.0416, +0.0517] | [+0.0412, +0.0520] | 1.07x |

**No conclusion changes: every headline still excludes zero.** That is exactly why the
correction is worth making. It costs nothing and removes an objection a referee would
otherwise land, and the gate `all_headlines_exclude_zero_clustered` means that if a future run
does push a headline across zero under clustering, the claim must be restated rather than the
narrower interval quietly retained. **The manuscript should print the protein-clustered
intervals throughout.**

## R1j: the expression confound, and why it cancels rather than explains

`expression_control.csv`, `expression_control_per_dataset.csv`, n = **84** of the 94 datasets
(the ten dropped have no usable expression track for their cell line).

**The objection.** Unbound genomic windows can be unbound for a boring reason: the gene is not
expressed in that cell line, so no RNA was there to crosslink. If the two arms' negative pools
differ in how untranscribed they are, the contrast could be an expression artifact rather than
a composition-matching one.

**The measurement, and the answer is balance.** **40.1%** of GC-arm negatives sit in a locus
with no gene above 1 TPM on either strand. The dinucleotide arm is at **40.1%** as well, and
the paired arm difference is **+0.0000122 [−0.0033, +0.0032]**, Wilcoxon **p = 0.64**. A
confound that is present in both arms in the same proportion cannot manufacture a difference
between them. That single fact is most of this section.

**But "cannot manufacture it" is not "does not touch it", so it was measured.** Restrict to
pairs whose negatives are transcribed and the contrast falls from this panel's **+0.0414** to
**+0.0310**, a drop of **−0.0104**. Nearly all of that is the restriction itself, not
expression: a **placebo stratified on region marginals**, dropping the same number of pairs at
random from the same region mix, loses **−0.0061** with no expression involved. The
pre-committed primary estimate is the difference,

> **expression-specific excess = −0.0043 [−0.0069, −0.0018]**

leaving an expression-corrected contrast of **+0.0371** and **89.7%** of this panel's contrast
surviving.

**Why the placebo is stratified, and it matters.** The unstratified placebo gives −0.0036;
stratifying on region moves it to −0.0043, so **+0.0007** of the naive estimate was locus mix
rather than expression. The restricted set is not a random subsample of the panel -- transcribed
negatives are enriched for exactly the regions where the contrast is largest -- and a placebo
that ignores that would charge expression for a compositional shift.

**One seed would have got this wrong, and did.** With five placebo seeds the interval covered
zero and the honest reading was "no expression artifact". With twenty it is −0.0043 with the
interval clear of zero. The per-arm Monte Carlo noise is **sd 0.0043**, the same size as the
effect, which is precisely why. This is the second time in this project that placebo-seed count
decided whether an artifact existed (R1c was the first), and it is a methods lesson worth one
sentence in the paper.

**Scope, and it must be stated.** The 89.7% is against **this panel's own** contrast of
**+0.0414** on 84 datasets, not against the published +0.0397 on 94. Comparing a subset estimate
to a full-panel published mean is the error that produced a retracted claim earlier in this
project, and the gate `r1j_expression_control` pins the panel's own contrast so the comparison
cannot silently drift.

**A second, larger number that is not the confound.** **67.6%** of GC-arm negatives are not
transcribed on the *labelled* strand, against 64.2% for the dinucleotide arm. That is the
strand criterion, not expression, and its effect is already measured and charged in R1c
(−0.0055, 85.4% surviving). The two controls overlap and must not be added.

## R1k: a third protocol, and the prediction it falsified

`three_arm_contrast.csv`, `three_arm_per_dataset.csv`, n = 94, 4-mer model only. Figure **f10**, which is now the paper's central figure.

Two protocols invite the reply that both are variants of one flawed design: composition-matched
unbound genomic windows, carrying everything R1c and R1j say about them. So a third was built
from data already on disk, implementing Horlacher et al. 2023's bias-aware **negative-2**:
negatives are **other RBPs' binding sites in the same cell line**, 1:1, excluding any near the
target's own sites. Those negatives are transcribed, CLIP-accessible and strand-correct by
construction, and no composition matcher touches them, so R1c, R1j and the sampler-hyperparameter
objection cannot apply to this arm at all.

| protocol | composition alone | apparent AUROC | **nested contribution** |
|---|---|---|---|
| dinucleotide-matched | 0.6274 | 0.6937 | **+0.0663** |
| GC-matched | 0.7827 | 0.8092 | **+0.0265** |
| **neg2** (other RBPs' sites) | **0.8248** | 0.8370 | **+0.0122** |

| contrast | value | positive in |
|---|---|---|
| dinuc − GC | **+0.0398** [+0.0337, +0.0459] | 88/94 |
| neg2 − GC | **−0.0143** [−0.0191, −0.0097] | 30/94 |
| neg2 − dinuc | **−0.0540** [−0.0634, −0.0447] | 6/94 |

**The prediction was written down before the run and it was wrong.** `three_arm_contrast.py`
predicted composition would fall toward 0.5 under neg2, because both classes are real crosslink
sites, and that the contribution would rise. Neither happened: composition is the **highest** of
the three arms under neg2, and the contribution the **lowest**. The reason is obvious in
hindsight and is a real finding: different RBPs bind compositionally different sites, so telling
one protein's sites from another's is largely a composition task.

**What this does to the paper.** It strengthens the thesis and kills the recommendation.

*Strengthens:* the same model, the same positives, the same folds and the same estimator give a
5.4-fold range in measured contribution. That is a much harder result to explain away than a
two-point contrast between two variants of one design.

*Kills:* "prefer dinucleotide matching, it reveals more of the model's contribution" is
**withdrawn**, and the reason is subtler than it first appeared.

**A correction I got wrong and a referee caught.** An earlier version of this section claimed
the ordering is NOT monotone in negative-set hardness, because the "bias-aware" protocol
reveals the least. That is false. Ordered by *measured* difficulty rather than by how principled
the protocol sounds, all three orderings are perfectly monotone:

| | dinuc | GC | neg2 |
|---|---|---|---|
| composition alone | 0.6274 | 0.7827 | 0.8248 |
| apparent AUROC | 0.6937 | 0.8092 | **0.8370** |
| nested contribution | **0.0663** | 0.0265 | **0.0122** |

**neg2 is the EASIEST of the three discriminations, not the hardest.** Telling one protein's
sites from another's is an easier task than telling them from matched genomic background, and
the contribution falls accordingly. So "harder negatives reveal more" is exactly true, and the
recommendation dies for a different reason: what "harder" buys is headroom, and headroom is set
by the protocol, so the advice reduces to "choose the protocol that leaves your baseline
weakest", which is a statement about the benchmark and not about the model.

**And it shows the decomposition cannot be rescued.** With three arms there are six transplant
residuals rather than four, and **their sign follows the direction of transplant**: carrying an
increment from a high-baseline arm to a low-baseline one gives a positive protocol effect
(neg2 → dinuc, **+0.0452**), the reverse gives a negative one (dinuc → neg2, **−0.0259**). That
is R1h's non-identification made visible across three protocols instead of two.

**The circularity critique, measured rather than argued.** Pooling all 282 arm-datasets, the
nested contribution tracks the composition baseline at **Spearman −0.600** (p = 6e-29). Most of
what the protocol does is decide how much room the baseline leaves. The paper should say this
outright rather than defend against it: it is the mechanism, not an objection to be deflected.

**Scope.** The 4-mer only. The deep models were never trained on the third arm, and comparing a
k-mer on one protocol against SpliceBERT on another would be worse than not running it.

## R1l: protocol and baseline are not separable, and that is the finding

`baseline_confounding.csv`, 282 protocol-dataset cells.

**The strongest objection left.** R1k concedes the contribution tracks the composition baseline
at Spearman −0.60. A referee will say: then "the protocol changes the measured contribution" is
just "the protocol changes the baseline, and the baseline changes the contribution", which is
trivially true and not worth a paper. This section is the answer, and it concedes half of it.

**Run naively, the test agrees with the referee.** Pooled over all 282 cells,
`gain ~ cubic(composition baseline)` gives R² = **0.3955**; adding protocol dummies gives
**0.4065**, F(2,276) = 2.55, **p = 0.080**. With dataset fixed effects, 0.8243 → 0.8301,
p = 0.045. Knowing which protocol produced a baseline adds almost nothing to knowing the
baseline.

**But that test cannot answer the question, and why not is the point.** The three protocols
barely overlap in baseline:

| protocol | median | 10th–90th percentile |
|---|---|---|
| dinucleotide-matched | 0.610 | 0.564–0.724 |
| GC-matched | 0.780 | 0.652–0.931 |
| neg2 | 0.815 | 0.719–0.937 |

The intersection of those ranges is **0.0056 AUROC wide and contains 3 of 282 cells**. There is
essentially no region where two protocols can be compared at the same baseline, so no analysis
on this data can separate them — and none could, because **the protocol IS the operation that
sets the baseline**. They are not two correlated variables. Asking "is it the protocol or the
baseline?" is malformed.

**What nevertheless survives, so the answer is not simply "it is all the ceiling".** If
compression were the whole story, transplanting a model's own d′ increment across baselines
would reproduce the other arm's gain. It does not:

| pair | observed difference | explained by compression |
|---|---|---|
| GC → dinuc | +0.0398 | **21%** |
| GC → neg2 | −0.0143 | **40%** |
| dinuc → neg2 | −0.0540 | **52%** |

And a single constant d′ increment applied to each cell's own baseline predicts all 282 gains at
only **R² = 0.183**, with residuals still ordered by protocol (**dinuc +0.0247**, GC −0.0038,
**neg2 −0.0143**). The arms differ by more than the arithmetic of the ceiling. **How much more
is not estimable**, and the paper says so rather than producing a number it cannot defend.

**What the reader should take from it.** Not "the protocol effect is X". The claim is that the
measured contribution is a joint property of model and protocol which no decomposition on
observational data can split, and the practical consequence follows directly: report the
composition baseline measured under the same protocol, because it is the only summary of what
the protocol did that a reader can act on, and never compare contributions measured under
different protocols, because there is no common scale on which to do so.

## R1m: no rescaling recovers a protocol-free quantity. This is the result that earns the title

`scale_sweep.csv`, n = 94, all three arms, 4-mer model. Figure **f11**. **The single most
important section after R1k.**

**The objection this closes, and it is the one a methods referee will reach for first.** AUROC
is compressive near 1, so of course a gain measured against a 0.82 baseline is smaller than the
same gain against a 0.63 baseline. Rescale to something that is not bounded -- d′, logit,
arcsine -- and the protocol dependence should evaporate. If it did, this paper would be a note
about a scale artifact.

It does not evaporate. Eight monotone transforms of the same 282 cells:

| coordinate | dinuc | GC | neg2 | **fold range** |
|---|---|---|---|---|
| raw AUROC gain | 0.0663 | 0.0265 | 0.0122 | **5.42** [4.43, 6.58] |
| Somers' D gain | 0.1325 | 0.0530 | 0.0245 | **5.42** [4.43, 6.58] |
| excess-normalised, g/(comp − 0.5) | 0.7649 | 0.1274 | 0.0421 | **18.18** [13.73, 23.96] |
| arcsine increment | 0.1426 | 0.0657 | 0.0357 | **3.99** [3.35, 4.80] |
| cloglog increment | 0.1851 | 0.0839 | 0.0470 | **3.94** [3.29, 4.80] |
| d′ increment (binormal) | 0.2645 | 0.1353 | 0.0833 | **3.18** [2.66, 3.89] |
| logit increment | 0.3089 | 0.1700 | 0.1134 | **2.72** [2.25, 3.39] |
| **headroom-normalised, g/(1 − comp)** | 0.1736 | 0.1167 | 0.0867 | **2.00** [1.67, 2.46] |

> **The floor is 2.00x [1.67, 2.46], and the interval excludes 1.**

**Read the table three ways.**

*Somers' D is the control.* It is an affine map of AUROC, so its range is identical to 5.42 by
construction. It is in the table to prove the sweep is measuring the transform and not noise:
any coordinate reporting something other than 5.42 there would be a bug.

*Rescaling helps, and the paper should say so.* Every unbounded link cuts the range, and the
ordering (arcsine, cloglog, d′, logit) is the ordering of how aggressively each stretches the
upper tail. So compression is **real and is part of the story** -- roughly a factor of two of
the 5.4 is scale. The claim is not that compression is absent; it is that compression is not
all of it.

*One transform makes it far worse.* Normalising by the *excess over chance*, g/(comp − 0.5),
inflates the range to **18x**. Dividing by a quantity that itself varies with the protocol and
can approach zero is not a stabilisation. That column is in the table because it is the
normalisation a reader is most likely to invent, and it is the worst thing they could do.

**And the winner is the paper's own recommendation.** The minimum over all eight is the
headroom-normalised coordinate g/(1 − comp) -- which is the two-number report, expressed as one
number. R1q then shows the same fact from the reader's side. That the recommended coordinate is
independently the least protocol-sensitive of eight candidates is the strongest form of the
paper's practical claim, and it was not designed that way: the sweep was run to try to kill the
result, and it selected the fix.

**THE SCOPE LIMIT WAS NOT A SCOPE LIMIT, IT WAS A HOLE. Retracted and replaced by R1u.** The
paragraph here used to say "eight transforms is not all transforms; a reparameterisation could
in principle reach 1.0x". It does not take a hypothetical. **g/(1 − comp) is not a coordinate:
it is the p = 1 member of the one-parameter family g/(1 − comp)^p, and this sweep stopped
exactly at the member that is the paper's own recommendation.** Extend the family and the range
collapses:

| p | 0 | 0.5 | **1.0** | 1.25 | **1.544** | 1.75 | 2.0 |
|---|---|---|---|---|---|---|---|
| fold range | 5.42 | 3.43 | **2.00** | 1.48 | **1.005** | 1.34 | 1.95 |

At **p = 1.544 the three protocols agree on the panel mean to 0.5%**. So "no monotone rescaling
recovers a protocol-free quantity" is **false as literally worded**, and the "eight standard
coordinates a benchmarker would use" hatch does not save it, because g/(1 − comp)^1.544 has
exactly as much derivation as g/(1 − comp) -- namely none. Both are argmins of the same sweep.
**A claim cannot present its argmin as a discovery and then dismiss a slightly larger family's
argmin as unprincipled.**

**What replaces it is a stronger claim, and R1u is where it lives:** the equalising exponent is
a property of the **benchmark**. Ours is 1.544; Horlacher's is **3.649**; and our exponent
leaves their benchmark at **2.340x** against the **2.381x** it started from, i.e. it buys
nothing there. **There is no *transportable* rescaling** -- any normalisation strong enough to
equalise one benchmark must be refitted on the next, which defeats the purpose. That is a
testable statement about transportability rather than an untestable one about a floor, and the
test uses data this project did not build.

**And the 2.00x inference itself survives against a corrected null.** The floor was being
compared against 1.0, which is too generous: max/min over three noisy means is bounded below by
1 and biased upward. With **equal true arm means and the real between-arm dataset pairing
preserved**, the null has median **1.064** and 95th percentile **1.175**. The observed 2.002 is
well outside it. Report the null, not the bare 1.0.

**So the honest reading of this section.** For the fixed coordinates a benchmarker would
actually reach for, the range floors at 2.00x and that is far outside chance. But protocol
independence is reachable by a *fitted* transform, and the fit does not transport -- which is
the real content of "there is no protocol-free measure".

## R1n: is it the protocol, or the baseline it leaves? Mostly the baseline, and the exception is named

`protocol_or_baseline.csv`, 282 protocol-dataset cells, 94 datasets. Figure **f12**. **This is
the paper's thesis.** R1l argued the question was malformed. A hostile referee found two places to ask it
anyway, and the answer is more interesting than either "it is all the baseline" or "the
protocol matters".

**1. Variance accounting.** Regress gain on a cubic in the composition baseline, pooled over
all 282 cells, then add protocol dummies:

| addition | incremental R² |
|---|---|
| protocol label, given the baseline | **1.0%** |
| baseline, given the protocol label | **11.0%** |

An order of magnitude, in the baseline's favour.

**2. The natural experiment, which is the sharpest design in the paper.** neg2 usually raises
the composition baseline relative to GC matching, but in **27 of 94** datasets it *lowers* it.
If the protocol label carried the information, neg2's deficit would persist there. It reverses:

| stratum | n | neg2 − GC gain | neg2 higher in |
|---|---|---|---|
| neg2 raises the baseline | 67 | **−0.0212** [−0.0269, −0.0157] | 14/67 |
| neg2 **lowers** the baseline | 27 | **+0.0028** [−0.0000, +0.0063] | 16/27 |

Within-dataset Spearman(Δbaseline, Δgain) = **−0.664, p = 3e-13**. Whichever protocol leaves
the lower baseline gets the higher contribution, **regardless of which protocol that is**.

**3. The dinucleotide-versus-GC contrast, matched on baseline.** comp_dn < comp_gc in **94 of
94** datasets, so within a dataset the arm and the baseline are perfectly rank-confounded by
construction and the comparison has to borrow across datasets. Matching each dinucleotide-arm
dataset to the GC-arm dataset with the nearest baseline (within 0.02 AUROC, n = 41), the
headline **+0.0398 falls to −0.0087 [−0.0265, +0.0122]**. The paper's own primary contrast does
not survive matching on the baseline. That is a concession, and printing it is what makes the
rest credible.

### The part that was withdrawn, and the mechanism that replaced it

An earlier version of this section claimed the protocol label carries *essentially no*
information. **That is withdrawn.** R1l's "common support is 0.0056 AUROC wide, 3 of 282 cells"
was the **three-way** intersection, and I ran the matched design only on the pair where it is
impossible. Pairwise:

| pair | common support width | cells | incremental R² of the protocol |
|---|---|---|---|
| GC vs dinuc | 0.073 | 33/188 | **0.05%** |
| **GC vs neg2** | **0.212** | **130/188** | **3.40%** |
| dinuc vs neg2 | 0.006 | 2/188 | 1.31% |

The GC-versus-neg2 comparison is not merely possible, it is easy. And run there, under two
independent designs, a protocol-specific residual **does** exist:

| design | neg2 − GC at matched baseline |
|---|---|
| nearest-baseline matching | **−0.0081** [−0.0130, −0.0036] |
| within-dataset intercept at zero baseline shift | **−0.0043** [−0.0077, −0.0012] |

Both exclude zero. The headline 1.0% averages the pair where the question is unanswerable with
the pair where it is answerable, and reporting only the average would hide the answer.

**The mechanism, and it makes the whole paper cohere.** The baseline gradient is a property of
**composition-matched** negatives and dies for other-RBPs'-sites negatives:

| arm | within-arm Spearman(baseline, gain) |
|---|---|
| GC-matched | **−0.545**, p < 0.001 |
| dinucleotide-matched | **−0.462**, p < 0.001 |
| neg2 (other RBPs' sites) | **−0.122**, p = 0.242, n.s. |

So there are two protocol *families*, not three protocols. Inside the composition-matched
family the baseline is the whole story (0.05%). Across families it is not (3.40%), and the
residual is a real property of what "other RBPs' sites" means as a control. **This is also why
R1p did not replicate the strong form**: Horlacher's negative-2 is the same family, so their
data shows the same residual. R1p is corroboration of the corrected claim, not a limitation
of it.

**The claim, in the form that survives.** The composition baseline is *most* of what determines
the measurable contribution -- 11.0% against 1.0%, reversing in 27 of 94 datasets, and
dissolving the headline contrast under matching. Within a protocol family it is essentially all
of it. Across families a residual of about **−0.005 to −0.008 AUROC** remains and is named
rather than absorbed.

## R1o: where the composition baseline stops is an analyst's choice, and it is doing half the work

`baseline_order.csv`, `baseline_order_per_dataset.csv`, n = **30** size-stratified datasets,
all three arms.

**The sharpest objection anyone made about what the quantity *is*.** A bag of 4-mer counts
carries no positional information, and the composition baseline stops at order 2. So "the
nested contribution of a sequence model over composition" is, definitionally, *order-3-and-4
composition beyond order-1-and-2 composition*, and where the baseline stops is an arbitrary
choice. `nested.py` defends order 2 as "exactly what a dinucleotide shuffle preserves", which
is a good reason to choose it and not a reason to believe the choice is immaterial. Nobody had
measured it.

Refitting the baseline at order 3 (mono + di + tri + entropy, 82 features, standardised, the
same `_oof_scores` fitter and the same DeLong comparison the paper uses everywhere):

| arm | composition, order 2 → 3 | over order 2 | over order 3 | removed |
|---|---|---|---|---|
| GC-matched | 0.7617 → 0.7837 | +0.0288 [+0.0188, +0.0403] | **+0.0058** [+0.0023, +0.0101] | **80%** |
| dinucleotide-matched | 0.6145 → 0.6666 | +0.0650 [+0.0506, +0.0800] | **+0.0105** [+0.0057, +0.0161] | **84%** |
| neg2 | 0.8005 → 0.8127 | +0.0122 [+0.0076, +0.0174] | **+0.0015** [+0.0004, +0.0027] | **88%** |
| **R1 contrast (dn − GC)** | | **+0.0361** [+0.0291, +0.0429] | **+0.0048** [+0.0025, +0.0071] | **87%** |

**The order-2 column is not a re-estimate, it is the published number.** The script's order-2
baseline *is* `nested.composition_features`, so its gain must equal `three_arm_per_dataset.csv`
cell by cell, and that is asserted per dataset rather than assumed: **max |difference| =
9.9e-17**. Without that anchor the order-3 column would be uninterpretable, and the first
version of this analysis proves the point -- see below.

**What it costs the paper, and it belongs in the Discussion, not a footnote.** **Four fifths of
the 4-mer's measured contribution, and 87% of the headline contrast, is one order of
composition.** The phrase "sequence model" has to be qualified hard: what a 4-mer adds over an
order-2 composition baseline is, almost entirely, order-3 composition -- **not motif
recognition and not anything positional**. Every gain remains positive with an interval clear
of zero, so something survives, but the honest magnitude over a trinucleotide baseline is
+0.0058 on the GC arm, not +0.0288.

**What it does not cost.** The protocol dependence, which is the paper's actual claim:

| | order-2 baseline | order-3 baseline |
|---|---|---|
| fold range across the three protocols | **5.34x** [3.95, 7.71] | **7.16x** [3.65, 22.65] |

**The range does not collapse. Whether it widens is not resolved, and the paper must not say it
does.** The change is **[−1.89, +16.43] with P(≤0) = 0.18**. The order-3 range's denominator is
the neg2 order-3 mean of +0.0015, whose own interval spans a factor of six, so the ratio is
barely identified. The claim that survives is the weaker one: **raising the baseline order does
not drive the protocol dependence toward 1**, which is the same division of labour as R1m -- a
defensible analyst choice changes the magnitude by a factor of five and does not remove the
protocol dependence.

*(An earlier version of this section said the range "widens", and the version before that said
it was "unchanged". Both overstated what 30 datasets can say, and both were printed from
NaN intervals as though point-identified. The removed fractions 80/84/88% likewise have
overlapping intervals, so **their ordering is not resolved either** -- the "neg2 loses the most,
GC the least" mechanism is withdrawn.)*

**The per-dataset sign, which a ratio of means cannot show.** This is the sharpest caveat in the
section:

| arm | order-3 gain positive in | median per-dataset ratio | top 3 datasets carry |
|---|---|---|---|
| GC-matched | 18/30 | +0.146 | 51% of the positive mass |
| dinucleotide | 23/30 | +0.119 | 38% |
| **neg2** | **17/30** | **−0.020** | **51%** |

**On the neg2 arm the typical dataset retains nothing** -- the median per-dataset ratio is
*negative* -- and its positive panel mean of +0.0015 is carried by three datasets. That mean is
the fold range's denominator. So "every gain remains positive with an interval clear of zero" is
true of the **arm means** and must not be read as a per-dataset statement.

**Scope, and the subsample has a blind spot the representativeness check does not cover.** n = 30,
size-stratified, because an order-3 fit is 82 features x 5 folds x 3 arms per dataset. On
order-2 quantities it is representative: its own three-protocol range is **5.34x** against the
full panel's **5.42x**, and its order-2 R1 contrast is +0.0361 against +0.0398 on 94.

**But the order-3 retained fraction is size-dependent** -- Spearman(pairs, retained) =
**+0.354** (GC, p = 0.055), **+0.627** (dinuc, p < 0.001), **+0.432** (neg2, p = 0.017), so
bigger datasets keep more -- and `head(30)` after `iloc[::3]` deterministically discards the
**two largest datasets** in the panel (32,384 and 19,533 pairs against the subsample's maximum
of 13,395). **The removal fractions are therefore biased upward.** R1r bounds that bias on the
full 94 for the GC and dinucleotide arms (78.3% against 79.9%, and 81.2% against 83.8%), but
**neg2 has no n = 94 counterpart and it is the fold range's denominator.**

**A correction, and it is the reason this section is trustworthy now.** The first version of
this analysis built its own composition block: all 4 mono + all 16 dinucleotide frequencies,
unstandardised, into sklearn's default `LogisticRegression`. Three departures from the paper's
baseline at once -- a singular design matrix, no entropy column, and an L2 penalty applied to
features on a 0.06 scale, which crushes the fit. It **underfit both baselines and therefore
overstated both gains**, by 1.22x on the GC arm and **2.31x on neg2**, and it underfit order 3
worst of all, reporting that only ~50% was removed when the true figure is 80--88%. It reported
a fold range of 2.67x that was not the paper's 5.42x and was never flagged as a different
quantity. It passed review because its order-2 R1 *contrast* landed at +0.0399 against the
published +0.0398 -- a coincidence of two offsetting errors. **A re-implemented baseline must be
asserted equal to the original, per cell, or it is a different experiment wearing the same
name.**

## R1r: an order-3 baseline absorbs the same amount from every model. That is nearly all of the 4-mer's contribution and a quarter of SpliceBERT's

`baseline_order_models.csv`, `baseline_order_models_per_dataset.csv`, n = 94, both arms, the
three model classes of R1g, on R1g's own committed per-window scores. Figure **f13**. No GPU:
the CNN and SpliceBERT scores are the ones R1g already used, so this costs nothing but CPU
time.

**Why this had to be run.** R1o's collapse is close to tautological for a bag of 4-mers, whose
only information beyond trinucleotide frequency *is* order 4. If it held for every model class,
the reading would be devastating and general: everything this literature measures as "what the
model adds over composition" is one order of composition. **It does not hold for every model**
-- but the reason is not the one the first version of this section gave, and the difference
matters.

GC arm, composition **0.7827 → 0.8023**; dinucleotide arm, **0.6274 → 0.6791**:

| model | GC: order 2 → order 3 | share | dinuc: order 2 → order 3 | share |
|---|---|---|---|---|
| 4-mer | +0.0265 → **+0.0058** | 21.7% [16.9, 26.2] | +0.0663 → **+0.0124** | 18.8% [14.7, 22.8] |
| CNN | +0.0330 → **+0.0200** | 60.5% [56.2, 64.6] | +0.0860 → **+0.0508** | 59.1% [54.7, 63.6] |
| SpliceBERT | +0.0890 → **+0.0668** | 75.0% [72.6, 77.2] | +0.1754 → **+0.1213** | 69.2% [66.6, 71.7] |

### The share is a ratio, and I read its denominator wrong

An earlier version of this section led with those shares and concluded "the collapse is a
property of the 4-mer, not of sequence models" -- i.e. that the k-mer is uniquely **fragile**.
**That reading is withdrawn.** The absolute amount the order-3 baseline absorbs is nearly
identical across model classes:

| model | GC: absorbed | dinuc: absorbed |
|---|---|---|
| 4-mer | +0.0208 [+0.0163, +0.0252] | +0.0538 [+0.0448, +0.0633] |
| CNN | +0.0130 [+0.0098, +0.0165] | +0.0352 [+0.0279, +0.0433] |
| SpliceBERT | +0.0222 [+0.0182, +0.0263] | +0.0540 [+0.0458, +0.0624] |

**SpliceBERT − 4-mer, absorbed: +0.0015 [+0.0003, +0.0027]** (GC) and **+0.0002 [−0.0022,
+0.0030]** (dinuc). On the dinucleotide arm they are indistinguishable; on the GC arm the
difference is real but two orders of magnitude smaller than the shares imply. The CNN actually
loses the **least** in absolute terms, which no capacity story predicts.

**But "a constant" was also wrong, and a statistician caught the second version too.** Two
things break it.

*The CNN absorbs materially less than either, on both arms.* cnn − kmer = **−0.0077** [−0.0098,
−0.0057] (GC) and **−0.0186** [−0.0226, −0.0148] (dinuc), intervals clear of zero. A "fixed
quantity from whatever model is being measured" does not survive its own middle row, and no
capacity story predicts the middle rung losing the least.

*And the paper's own compression correction was never applied to the order step.* An order-3
baseline **raises** composition (+0.0196 GC, +0.0517 dinuc), cutting the headroom every model is
working in -- and it cuts **more** from a model with a bigger gain, which is exactly what made
the raw absorption look flat. `deep_model_contrast.py` calls this correction "not optional"
elsewhere in the paper. Transplanting each model's order-2 d′ increment onto the order-3
baseline:

| GC arm | predicted by compression alone | observed | **residual** |
|---|---|---|---|
| 4-mer | +0.0239 | +0.0058 | **−0.0182** |
| CNN | +0.0298 | +0.0200 | **−0.0099** |
| SpliceBERT | +0.0809 | +0.0668 | **−0.0141** |

**Corrected for the moving ceiling, the 4-mer does lose more -- by 1.29x (GC) and 1.41x
(dinuc), not by the 3.4x the shares imply.** So both of this section's earlier framings were
partial: the shares overstate the difference, the raw absorption understates it, and the honest
figure is about **1.3 to 1.4x**. The CNN still loses the least on the corrected scale, which is
gated.

**The statement that survives.** Raising the baseline one order absorbs **roughly +0.021 AUROC
on the GC arm and +0.054 on the dinucleotide arm**, only weakly dependent on model class; after
correcting for the headroom it removes, a bag of 4-mers loses about **1.4x** what a fine-tuned
transformer does. That is nearly all of the 4-mer's contribution and about a quarter of
SpliceBERT's -- but because the totals differ, not because the 4-mer is three times as fragile.

### What the result actually is, and it is still strong

*1. Over a trinucleotide baseline the 4-mer's contribution is not reliably positive; SpliceBERT's
is.* This is the per-dataset fact the shares cannot show, and it is the sharpest thing here:

| order-3 gain positive in | GC arm | dinuc arm |
|---|---|---|
| 4-mer | **65/94** | 75/94 |
| CNN | 80/94 | 90/94 |
| SpliceBERT | **94/94** | **94/94** |

**For 29 of 94 datasets a 4-mer adds nothing at all over trinucleotide composition** -- the
point estimate is negative. SpliceBERT is positive on every dataset in both arms. A mean of
+0.0058 was hiding a third of the panel at or below zero.

*2. The order-3 baseline is a sharper instrument than the order-2 one, and this is the useful
recommendation.* SpliceBERT's advantage over the 4-mer widens from **3.36x to 11.61x** on the GC
arm and **2.65x to 9.75x** on the dinucleotide arm, and at order 3 the gaps are clean:
SpliceBERT − 4-mer = **+0.0611 [+0.0509, +0.0720]** (GC) and **+0.1089 [+0.0984, +0.1194]**
(dinuc), CNN − 4-mer = **+0.0142 [+0.0104, +0.0184]** and **+0.0384 [+0.0311, +0.0468]**. A
benchmark that wants to distinguish model classes rather than rank them by how much
short-range composition they absorb should raise the composition baseline, not the model.

**And the recommendation must be stated with the circularity it contains, or it is rigged.**
Trinucleotide counts are **linear aggregates of 4-mer counts**, so an order-3 block absorbs a
4-mer's score partly *by algebra*. Measured as the R² of each model's standardised
out-of-fold score on the composition design, the increment from adding the trinucleotide block
is **+0.152** for the 4-mer against **+0.051** (CNN) and **+0.068** (SpliceBERT) on the GC arm.
**The order-3 ranking is close to the ranking of how orthogonal each model's score is to the
baseline's span** -- which is not the same thing as capacity. So the honest form of this
recommendation is: *raising the baseline's order separates models by how much of their signal
lies outside short-range composition, and it is structurally unfavourable to bag-of-k-mers
models because their features are nested inside the baseline's.* Stated that way it is useful.
Stated as a capacity instrument it is a circularity in method's clothing.

Both are corollaries of the constant: a fixed absolute absorption leaves a model with a large
contribution looking far better than one with a small contribution, and that IS what a
benchmark should want to distinguish. But the mechanism must be stated as absorption, not as
robustness.

*3. The paper's own claim survives at the higher baseline, for every model class.* The
protocol contrast at an order-3 baseline:

| model | order-2 contrast | order-3 contrast |
|---|---|---|
| 4-mer | +0.0398 [+0.0324, +0.0478] | **+0.0067** [+0.0044, +0.0096] |
| CNN | +0.0530 [+0.0438, +0.0630] | **+0.0309** [+0.0246, +0.0374] |
| SpliceBERT | +0.0864 [+0.0777, +0.0954] | **+0.0545** [+0.0488, +0.0602] |

All six exclude zero. **Protocol dependence is not an artefact of the baseline's order for any
model class**, and for the two models that survive the higher baseline it is barely reduced.

**And the withdrawn claim stays withdrawn.** On the ratio scale, dinuc/GC as a ratio of panel
means over 94 datasets: **2.50x / 2.61x / 1.97x** at order 2 and **2.16x / 2.55x / 1.82x** at
order 3. SpliceBERT is the **lowest** at both orders. The order-3 table must not be used to
revive "the contrast grows with capacity". *(These are a ratio of panel means on n = 94 and are
deliberately a different estimator from R1g's 3.08 / 3.51 / 2.38, which is a geometric mean of
per-dataset ratios on the 77 datasets positive in both arms for all three models. Both are
reported because neither is wrong; do not present them as the same quantity.)*

**Scope and cost.** No new training. These are R1g's committed per-window CNN and SpliceBERT
scores, on R1g's intersected row set, imported from `deep_model_contrast.py` rather than
reimplemented. The order-2 column is therefore R1g's published column, asserted per cell at
**max |difference| = 1.0e-16**. Two arms only: the deep models were never trained on neg2.

## R1s: whose property is the multiplier? The protein's, but the comparison has to be fair

`multiplier_variance.csv`, 262 (dataset x model) cells with both arms positive, 79 proteins.

**Why this section exists.** The README carried "the multiplier is a property of the PROTEIN
(68.9% of variance) not the model (1.5%)" as R1g's headline. That decomposition was run during
review and **committed to no script and no table** -- the most quotable line about R1g was
sourced to a conversation. It is now computed, gated, and corrected in one respect that matters.

Decomposing the **log** multiplier log(gain_dn / gain_gc), each factor's share taken as its
drop in residual sum of squares given the others:

| factor | levels | share | **null, relabelled data** | **excess** | p |
|---|---|---|---|---|---|
| **protein** | 79 | **64.8%** | **29.7%** | **+35.1** | < 0.0005 |
| model class | 3 | 2.8% | 0.8% | +2.0 | **0.023** |
| cell line | 2 | 0.2% | 0.4% | −0.2 | 0.49 |
| residual | | 32.8% | | | |

Panel multiplier **2.99x**, as exp of the mean log. *(A ratio of panel means is a different
number and R1g reports a third, a geometric mean over the 77 datasets positive in both arms
for all three models. Name the estimator every time.)*

**The first correction, and it is a fairness correction.** "68.9% against 1.5%" is not a
comparison between factors. **A factor with 79 levels absorbs 29.7% of the variance of
*relabelled* data**, while a 3-level factor absorbs 0.8%, so most of the apparent gap is
degrees of freedom. The comparable statistic is each factor's excess over its own null.

**The second correction, found by an adversarial statistician, and it cost four fifths of the
headline.** That wholesale-permutation null is **the wrong null**. Permuting *y* destroys the
protein effect *and* the within-dataset correlation among the two or three model cells that
share the same windows, labels and folds. And because **79 proteins sit on 94 datasets, protein
is very nearly the *dataset* factor** -- entered on its own, dataset explains **68.0%**, *more*
than protein's 64.8%. So the wholesale null is far too permissive:

| null | null share | observed | **excess** | p |
|---|---|---|---|---|
| wholesale *y* (as first published) | 29.7% | 64.8% | **+35.1** | < 1/2001 |
| **protein labels permuted between datasets, blocks intact** | **57.6%** [49.8, 63.6] | 64.8% | **+7.2** | **0.0030** |

**The claim survives and its magnitude does not.** The honest reading: 64.8% is a
*dataset*-level share, and the part attributable to **protein identity rather than to the
individual experiment is +7.2 points**, p = 0.0030.

**And the section now leads with a direct test instead of a share**, because it is stronger and
was sitting unused: does the *same protein* get the same multiplier in the *other cell line*?

> **r = +0.586, Spearman +0.594 (p = 0.0001), over 40 protein x model pairs** from the 15
> proteins assayed in both K562 and HepG2. An ICC of about 0.59.

That is the evidence for "a property of the protein". The variance share is a weaker restatement
of it.

**Two further corrections to the wording.**

1. **Cell line is not "noise", it is *not detectable*.** Only 15 of 79 proteins appear in both
   lines, so the conditional cell effect is identified off 15 proteins. p = 0.49 is a valid test
   of a hypothesis this design has almost no power to reject. Say "not detectable, with 15
   informative proteins".
2. **Model class is not null** (p = 0.023), and this was checked against the censoring visible
   in its own input: "both arms positive" drops 14 CNN cells, 6 k-mer and **0 SpliceBERT**, so
   the model comparison is censored *by model*. On the balanced 77-dataset panel the model share
   *rises* to 4.1% (p = 0.008) and the ordering is unchanged. The claim is robust; say that it
   was checked.

**A caveat on the estimand itself.** The log multiplier is mechanically dominated by its
denominator: Spearman(gain_gc, log multiplier) = **−0.720**, and five quintiles of the GC-arm
gain alone absorb 38.4% of the variance given model and cell. So part of what looks like a
"protein property" is the protein's effect-size *scale* plus small-denominator inflation, not
protocol biology. State it.

**The gates.** `null_must_be_reported` fails the build if the decomposition is published without
a null; `block_null_must_be_reported` additionally requires the **smaller** excess to be present,
so the five-fold-too-generous figure cannot be quoted on its own again.

## R1p: external validation, and where it limits us

`horlacher_arm.csv`, n = 45 of the 88 datasets shared with Horlacher et al. 2023's release.
Zenodo `10.5281/zenodo.10600977`, md5 verified. **Their positives, their peak calling, their
negative sets, their fold assignments.** Only the measurement is ours.

| their negative set | composition | apparent | **contribution** |
|---|---|---|---|
| negative-1 (transcript background) | 0.8211 | 0.8606 | **+0.0395** |
| negative-2 (other RBPs' sites) | 0.7560 | 0.7726 | **+0.0166** |

**What replicates, and it is the strongest external support in the paper.** The fold range on
their benchmark is **2.38x**, against **2.50x** for our corresponding two-protocol comparison.
Measured on data this project did not build, with folds it did not choose. And the baseline
gradient replicates in sign: within-dataset Spearman(Δbaseline, Δgain) = **−0.372, p = 0.012**,
against our −0.664.

**What does not replicate, and it limits R1n.** On our data the natural experiment reverses:
where the alternative protocol *lowers* the baseline, the contrast flips (−0.0212 → **+0.0028**).
On theirs it does not. The deficit shrinks in the predicted direction — **−0.0409** where
negative-2 raises the baseline, **−0.0184** where it lowers it — but stays negative in both
strata. Their negative-2 gives less contribution regardless of which way the baseline moved.

**So the honest claim is the weaker one.** "The protocol label carries essentially no
information beyond the baseline" holds on our windows and **not** on an independent benchmark.
What survives everywhere is that the baseline explains most of the protocol dependence, and
that the range itself travels. Both the replication and the failure are gated, so a future run
cannot quietly restore the strong form.

### The stronger result that was sitting in this table unreported

This section described the external data as a partial *failure* that limits R1n. **That
understated it.** R1n's real content is the **two-family mechanism**: the baseline gradient is a
property of composition-matched negatives and absent for other-RBPs'-sites negatives. Tested
per arm on their data (R1u):

| their arm | family | within-arm Spearman(baseline, gain) |
|---|---|---|
| negative-1 (transcript background) | composition-matched | **−0.645, p = 1.8e-06** |
| negative-2 (other RBPs' sites) | other-RBPs'-sites | **−0.187, p = 0.22, n.s.** |

Against ours: **−0.545 / −0.462** for the two composition-matched arms and **−0.122 n.s.** for
neg2. **The mechanism replicates, arm for arm and family for family, on an independent lab's
negatives and fold assignments, at p = 2e-06.** This is the strongest external result in the
paper and it takes five lines against a table that was already committed.

It also **explains** the arm-level inversion this section previously reported without comment.
Their negative-2 has the lower baseline *and* the lower gain, which is anti-monotone to our
three-arm ordering -- and that is exactly what a cross-family residual predicts. So R1p is
**corroboration of the corrected R1n, not a limitation of it.**

**And the range comparison should be family-matched.** Their negative-2 *is* our neg2, so the
corresponding pair in our data is **GC vs neg2 = 2.17x**, not the dinuc-vs-GC **2.50x**
previously quoted against their 2.38x. The family-matched comparison is both more defensible
and closer.

## R1q: the recommendation is tested, not asserted

`recommendation_works.csv`, n = 94, intervals bootstrap the 79 proteins.

Everything above shows that **not** reporting the composition baseline is a problem. Nothing
showed that **doing** it helps, and a benchmarking paper whose recommendation is untested is a
complaint rather than a contribution. So: put a reader in the position the recommendation is
meant to rescue -- two papers, two protocols, comparing what the models contributed -- and
compare raw contributions against the headroom-normalised coordinate the two-number report buys.

| protocol pair | rank agreement, raw → headroom | disagreement on a common scale |
|---|---|---|
| GC vs dinuc | 0.848 → **0.899** [+0.007, +0.115] | 1.151 → **0.482** |
| GC vs neg2 | 0.674 → **0.724** | 0.504 → **0.451** |
| dinuc vs neg2 | 0.607 → **0.628** | 1.552 → **0.802** |

**Rank agreement improves in 3/3 pairs and disagreement shrinks in 3/3**, by 58%, 10% and 48%.
Rank agreement is scale-free, so no normalisation can flatter it; the disagreement column
divides each coordinate by its own panel mean, because raw and headroom-normalised numbers have
different units and comparing their absolute differences directly would rig the result.

**The limit, gated so the wording cannot drift.** Only the GC-versus-dinucleotide rank
improvement has an interval clear of zero. The other two are consistent in direction and not
individually significant. The claim is therefore **direction-consistent**, not uniformly
significant, and the disagreement reduction is the stronger half of it.

**What would have falsified it.** Rank agreement falling, or disagreement rising, under
normalisation -- in which case the honest paper says "we can show the problem and cannot offer
a fix", which is still publishable and considerably weaker.

### It falsifies out of sample. Both criteria fire, and the paper must say so

Run the identical test on Horlacher's 45 datasets -- their positives, negatives and folds
(R1u, `transport_check.csv`):

| | rank agreement | scale-free disagreement |
|---|---|---|
| raw | **+0.706** | **0.860** |
| headroom-normalised | **+0.656** | **0.908** |

**Rank agreement falls and disagreement rises.** Those are, verbatim, the two criteria
`recommendation_works.py` names as falsifying. The change in rank agreement is **−0.050
[−0.222, +0.140], P(≤0) = 0.70**, so at n = 45 this is a **failure to replicate rather than a
refutation** -- but it is the paper's own pre-registered test applied to the only data the
project did not build, and it points the wrong way.

**And the in-sample result is weaker than it looked, for a structural reason.** The only pair
whose rank improvement clears zero is **GC versus dinucleotide** -- which is exactly the pair
R1n shows is **0.05% protocol and essentially all baseline**. The fix is validated where it
cannot fail and is null in the two pairs that cross protocol families. The mechanism is visible
in R1p's own table: their negative-2 has the **lower** baseline (0.7560 vs 0.8211) *and* the
lower gain (0.0166 vs 0.0395), so at the arm level the headroom relation **inverts** there.

**The honest claim, and it is a downgrade.** *We can demonstrate the problem; we cannot offer a
normalisation that is shown to transport.* Reporting the composition-only AUROC under the same
protocol remains the right practice -- it is what makes the problem visible and it is what the
two-number report is for -- but the paper must **not** claim that normalising makes
cross-protocol contributions comparable. R1l already says they are not comparable; the strong
form of R1q contradicted R1l, and R1l wins.

## FOLD INTEGRITY: 20 of the 94 dinucleotide-arm deep score sets did not use the study's folds

`fold_integrity.csv`, `fold_integrity_per_dataset.csv`, 376 (dataset, arm, model) score sets.
**The most serious defect found in review, and the harness could not see it.**

**What happened.** `config/folds.tsv` is a frozen chromosome-to-fold map, and chromosome
grouping is the entire design: a model must never be tested on a chromosome it trained on. For
**20 of the 94 dinucleotide-arm datasets**, the committed CNN and SpliceBERT per-window scores
were produced under a *different* partition -- a stratified random split that **preserved fold
sizes**, so no count or balance check could detect it, while putting up to **23 chromosomes in
every fold**.

**Confirmed three independent ways.** The score file's own `fold` column disagrees with
`dataset.tsv`'s (agreement 0.48 to 0.99); every affected score fold spans all 23 chromosomes
where every table fold spans 4 or 5; and the direct metric -- the fraction of rows having a
same-strand genomic neighbour within 1 kb assigned to a **different** fold -- reaches **44.5%**
(CSTF2T:HepG2), against **exactly 0** under `folds.tsv` and **0 everywhere in the GC arm**,
because under chromosome grouping a near neighbour is necessarily in the same fold.

| arm | chromosome-grouped | aligned to folds.tsv | max chromosomes in a fold | max cross-fold 1 kb neighbours |
|---|---|---|---|---|
| **GC-matched** | **94/94** | **94/94** | 5 | **0** |
| dinucleotide-matched | **74/94** | 74/94 | **23** | **44.5%** |

**How it survived 549 assertions.** `deep_model_contrast.py` asserted *in its own docstring*
that "both arms use the same folds, the same seed... the only difference is how the negative
windows were chosen." That sentence was false for those 20 datasets and nothing tested it. The
harness gated the **values** the scores produced and never **what produced them** -- which is
this project's oldest and most expensive lesson, now recurring at the level of the evidence
files rather than the summary tables.

**What moves, and it is less than the defect suggests.** The k-mer is refit in-script on the
study folds for every dataset, so it is a clean internal control. Dropping the 20:

| R1 contrast (dn − GC) | all 94 | **74 chromosome-grouped** | shift |
|---|---|---|---|
| 4-mer | +0.0398 | **+0.0436** | +0.0038 |
| CNN | +0.0530 | **+0.0506** | −0.0024 |
| SpliceBERT | +0.0864 | **+0.0878** | +0.0015 |

**Every shift is inside the protein-clustered half-widths of about 0.009**, so this is a
disclosed limitation and not a correction to any claim. The quantity that *is* affected is the
**CNN's dinucleotide-arm deep gain**: a difference-in-differences against the k-mer gives
**+0.0294** on the leaky datasets, though they are also the panel's largest (median 15,596
against 4,373 pairs) so that figure is an upper bound confounded with size.

**And the scope of what is untouched is the important part. The GC arm is clean, 94 of 94, and
every k-mer number in this paper is refit on the study folds.** So R1's headline, R1m, R1n,
R1q, and the GC arm of R1o and R1r are unaffected. What must carry the caveat: R1g's
dinucleotide-arm deep gains, R1r's dinucleotide arm, and R1s, whose panel multiplier moves from
2.99x to about 2.83x when the affected datasets are dropped.

**The gates now enforce three things**: the GC arm must stay perfectly grouped with a
cross-fold neighbour fraction of exactly zero; the leaky count is pinned two-sided at 20, so a
rerun that leaked *more* cannot pass; and every contrast must survive dropping them, or the
dinucleotide deep sweep has to be rerun rather than caveated.

## R1t: the ladder's rungs were not measured on the same scale. Measured, and the claim survives

`score_scale_check.csv`, `score_scale_per_dataset.csv`, n = 94, both arms.

**The defect, found by a methods referee and confirmed in full.** The nested contribution adds
the model's score as a covariate to a logistic regression. The 4-mer entered as a **log-odds**
(`baseline.oof_scores` returns `decision_function`); the CNN and SpliceBERT entered as
**probabilities**, because that is what their GPU sweeps wrote to `data/evidence`. A logistic
regression is **not** invariant to a monotone *nonlinear* transform of a covariate, so the deep
models' gains were measured with a range-compressed version of their own score while the
k-mer's were not. AUROC itself *is* invariant to monotone transforms, so no standalone AUROC in
this paper is affected -- but the nested fit can move, and this is an asymmetry inside the one
section whose entire point is that the three rungs share rows, folds and estimator.

**Measured on both scales, all 94 datasets.** The probability-scale refit reproduces the
published deep gains at **9.97e-17**, so what follows is the scale and not a different pipeline.

| | as published | on logit(p) | **difference** |
|---|---|---|---|
| CNN, GC | +0.0330 | +0.0338 | **+0.0008** [+0.0006, +0.0011] |
| SpliceBERT, GC | +0.0890 | +0.0920 | **+0.0030** [+0.0023, +0.0037] |
| CNN, dinuc | +0.0860 | +0.0866 | **+0.0006** [+0.0004, +0.0009] |
| SpliceBERT, dinuc | +0.1754 | +0.1787 | **+0.0033** [+0.0024, +0.0042] |

**The bias is real, one-directional and small.** All four intervals exclude zero, so this is a
systematic effect and not noise. Its direction is that **the probability scale understates the
deep models**, which makes every published deep number in this paper *conservative*. The
largest correction is +0.0033, about 2% of SpliceBERT's dinucleotide-arm gain.

**And the paper's claim is invariant to it, which is what makes this a disclosed sensitivity
rather than a correction:**

| R1 contrast (dn − GC) | as published | all-logit |
|---|---|---|
| 4-mer | +0.0398 | +0.0398 (unchanged, it was always the logit) |
| CNN | +0.0530 | **+0.0528** |
| SpliceBERT | +0.0864 | **+0.0867** |

Shifts of **0.0002 and 0.0004**, against protein-clustered half-widths of ~0.009. The ladder's
order, the contrast and the multiplier all stand.

**What the manuscript does.** Report the published values, state the asymmetry, and give the
measured correction, gated by `max_contrast_shift` so no future run can quietly let it grow.
The alternative -- silently refitting the deep arms on the logit and restating gated numbers
that appear throughout this document -- is precisely the move this project has caught itself
making before. **Measure first, then decide; and say which one was published.**

**Two related asymmetries I have NOT closed, and they belong in Limitations.** First, the
composition features and the score are standardised over **all** rows including the held-out
fold before the out-of-fold fits: unsupervised and small, but it is a use of test-fold
statistics in a paper whose discipline is "out-of-fold on both sides". Second, the nested full
model is an out-of-fold **stack**: the score added for a training row in fold *g* came from a
model trained on folds excluding *g*, which includes the test fold *f*. The composition-only
arm has no fitted feature and so no such leak, making the asymmetry point in the direction
that **inflates** the nested contribution. R1i's design effect of 1.35 addresses the *variance*
consequence of fitted scores; nothing here addresses the *bias*. Both are stated rather than
bounded.

## Match quality: what the matchers achieved, not what they were asked for

`match_quality.csv`, `match_quality_per_dataset.csv`, all **456,734** pairs per arm.

**Why this section exists.** Every arm is specified by its tolerance or its objective, and this
paper argues *from* the specification: the GC arm controls **1 of 15** degrees of freedom of the
composition baseline and the dinucleotide arm controls **15 of 15**, which is why the sign of the
headline contrast is design-implied. That is exactly true of the *design*. The *realisation* was
measured nowhere -- no committed table, no golden key, and none of the previous 522 assertions
touched it. `prepare.py` did compute the L1, into `prep_report.txt` files that no longer exist.

| arm | median \|ΔGC\| | p90 | p99 | max | median dinuc L1 | p90 |
|---|---|---|---|---|---|---|
| GC-matched | **0.0297** | 0.0495 | 0.1089 | **0.1486** | 0.500 | 0.780 |
| dinucleotide-matched | 0.0198 | 0.0693 | 0.1485 | 0.4159 | **0.220** | 0.340 |
| neg2 | **0.1387** | 0.3169 | 0.4555 | 0.6733 | 0.700 | 1.060 |

**Two facts the Methods must state, and neither was being stated.**

*1. The GC tolerance is not 0.05.* The matcher relaxes: 40 draws, tolerance **doubled to 0.10**
after 25 misses, then a best-seen candidate accepted up to **3x = 0.15**. **94.8%** of pairs are
inside the nominal 0.05, 98.7% inside 0.10, 100% inside 0.15 -- and the observed maximum,
**0.1486**, sits exactly at the fallback cap. Quoting `gc_tolerance: 0.05` overstates the match
for one pair in twenty. Report the ladder.

*2. "Dinucleotide-matched" means L1 ≈ 0.22, not ≈ 0.* The matcher improves median composition
distance from **0.500 to 0.220**, a factor of **2.27**, in frequency units on a 0-2 scale. So
"15 of 15 degrees of freedom controlled" is exactly true of the design and **approximately**
true of the data. This is also the number R1f needed: that section is described as "confounded
with achieved match quality" and until now no figure for it existed anywhere.

**And neg2 is not composition-matched at all** -- median \|ΔGC\| **0.1387**, 4.7x the GC arm,
with only 20.4% of pairs inside 0.05. That is by construction, and it is precisely why its
composition baseline is the highest of the three arms (0.8248). The three protocols are not
three tolerances of one matcher; they are three different operations, and this table is the
evidence for that framing.

## Supplementary Table S1: the panel, joined to ENCODE

`supplementary_table_s1.csv`, built by `scripts/table_s1.py`, wired into run.sh stage 13b.

One row per dataset: protein, cell line, **ENCFF file and ENCSR experiment accession**,
replicate count, pair count, whether the dataset carries all three protocols, and each
protocol's composition baseline and nested contribution. **95 datasets, 79 proteins, 95
distinct ENCODE experiments, 2 replicates each.**

**It also settles the 94-versus-95 question**, which appears in several tables and must be
stated once in Methods rather than explained per table: **94 datasets carry all three
protocols and are the unit of every three-arm result; 1 is ladder-only** (R2, n = 95). The
script fails loudly on an incomplete accession join rather than emitting a table with silent
gaps, because a supplementary table with holes looks complete.

## R2: the model ladder (methods table, not a result)

`matched_four_models.csv`, n = 95, identical chromosome-level folds. Figure **f2**.

composition **0.6279** < k-mer **0.6875** < CNN **0.7063** < SpliceBERT **0.8091**.
SpliceBERT beats composition on **95/95**.

**Prior art:** Horlacher et al. 2023, *Briefings in Bioinformatics* published the negative-set
effect across 11 RBP methods. This is scaffolding for R1, cited as replication.

## Not in this paper, and why

Six sections that appeared in earlier drafts are cut. They are listed because a reader of the
repository will find their code and tables still present, gated and reproducible, and should
know they were considered and rejected rather than overlooked.

| section | why it is gone |
|---|---|
| **R2**, the four-model ladder | Horlacher et al. 2023, *Brief Bioinform* 24(5):bbad307 published the negative-set effect across 11 methods and 313 experiments. Retained above only as a methods table and cited as replication. |
| **R3**, pooled vs paired AUROC | Prior art (van Klaveren et al. 2016, *Stat Med*; Janes & Pepe 2008), and on inspection a third view of the composition story rather than an independent result. |
| **R4/R4b/R4c**, the ClinVar ladder | Retracted. Every fit passed conservation into the model, so the row labelled "controls = none" was already adjusted and the published "0.21% attenuation" measured the removal of 168 rows. The AUROC framing it replaced is separately circular: ClinVar pathogenic assertions lean on PP3 evidence from conservation-derived tools, so ranking a model against phyloP partly ranks it against the answer key. |
| **R4d**, the corrected attenuation analysis | See the limitation below. The remedy is Schuster et al. 2021, *BMC Med Res Methodol* 21:136, and with R4 cut there is nothing left in this paper for it to correct. |
| **R5**, the wrong-protein specificity control | Conditional and underpowered. The detection plateau moved from 20 to at least 45 pathogenic variants under scrutiny, and the downsampling experiment that would settle it was never run. |
| **R6**, the substitution-spectrum baseline | Its confidence interval is wider than half the effect it reports, and it spans zero. |
| the composition-**share** framing | An algebraic identity: share_m = C/gain_m with C constant across models, so the model contrast excluded zero with probability 1. |

**One limitation inherited from the cut ClinVar work, stated because it was paid for.** A
conservation control that compares a coefficient before and after adjustment must be read
against a calibrated null and not against zero: logistic regression is not collapsible, so a
strong covariate that is genuinely independent of the exposure *amplifies* the other
coefficient rather than leaving it alone. An analysis in this repository originally read a
near-zero attenuation as evidence of independence when it is evidence of the opposite. The
argument is general and belongs to Schuster et al. 2021; the specifics are in
`scripts/unconditional_refit.py` and are not claimed here.

## RETRACTED, on the record

**The "model-dependent composition share."** Composition reproduces 0.682 of a k-mer model's
gain, 0.620 of a CNN's, 0.414 of SpliceBERT's, contrast +0.268. It is an **algebraic identity**:
`share_m = C / gain_m` where `C = mean(composition_auroc) − 0.5 = 0.127853` is identical across
models, so `share_kmer / share_SB == gain_SB / gain_kmer` exactly, **1.648166 both ways to 6
decimal places**. It is a monotone rescaling of R2, which is prior art, and SpliceBERT beats the
k-mer model on 95/95 datasets, so the contrast excluded zero with probability 1. `verify.py`
now gates the identity rather than a confidence interval. It appears in the manuscript as one
footnote warning against share-as-comparator estimators.

**The pooled ClinVar ladder (matched 0.829).** Simpson-inflated: on a fixed 82-dataset panel
the pooling inflation is **+0.165**, against +0.032 for conservation and +0.021 for the k-mer
model, the deep arm is the one pooling flatters most. Supplement only, labelled as such.

**R3-locality (ISM positional concentration).** Cut: architecture-confounded by construction.

**The "clean donor" stratification (+0.121, 16/17).** The stratifier `shared_frac = |t∩d|/|t|`
is normalised by the target, correlates with donor variant-set size at rho **+0.696**, and the
flagship association reverses and dies under adjustment (+0.299 p=0.049 → −0.240 p=0.12).

---

## LIMITATIONS: manuscript-ready

**1. Negative strand assignment.** Sampled negatives inherit the *positive's* strand while
their location is chosen to match composition, so strand is close to a coin flip. Audited on 40
datasets: **55.2%** of negatives carry the strand of the gene they sit in among *unambiguous*
windows; a further **14.0%** lie in genes transcribed on both strands, so the demonstrably-sense
fraction is **47.4%**. None are intergenic. Roughly half the negative set is therefore antisense
sequence no transcript produces, while every positive is true sense RNA. **This inflates all
absolute AUROCs reported here**, in both arms, and no absolute number in this paper should be
compared against a published benchmark that samples negatives differently.

R1c measures what it does to the *contrast*, which is the quantity the paper claims, and bounds
it at **−0.0055** [−0.0089, −0.0022] with **85.4%** surviving. Two things limit that. The
placebo design is an assumption, not a proof: matching on region and gene density removes a
locus mix that turns out to be **+0.0004** [−0.0010, +0.0019], an interval straddling zero, so
that component is seed noise rather than a measured confound, and the design cannot remove a
confound we did not think to match on. And the estimate rests on 40 of the 94 datasets, those
with a strand audit and canonical window tables in both arms.

*Correction, found by a referee and not by the harness.* This paragraph previously read
**−0.0036** [−0.0071, +0.0001] with **90.6%** surviving, and quoted a locus-mix component of
**−0.0024** with an "interval clear of zero". All three are the 5-seed estimates that the R1c
section above explicitly retracts: at 20 placebo seeds the excess is −0.0055 and the locus mix
collapses to noise. The stale values sat in the section labelled "manuscript-ready" while the
correct ones sat 280 lines earlier in the same file, and **`verify.py` passed 284/284 and
`audit_manuscript.py` reported no orphan**, because a retracted number is still traceable to a
committed table -- the older `strand_placebo.csv` rows it was computed from. Traceability is not
currency. See the lesson added to `docs/61`.

*An earlier draft defended this with the composition-share contrast, which is retracted below as
an algebraic identity, and one of the numbers it quoted appeared in no committed table. That
defence is withdrawn and R1c replaces it. A further sentence claiming this strand convention is
shared with the matched-negative samplers this work benchmarks against is also deleted: we have
no citation for it and it is probably false, since samplers such as GraphProt draw unbound
windows from annotated transcripts in transcript space rather than genome space.*

**1a. The sign of the contrast is implied by the design; only its magnitude is a result.**
Stated in R1 and repeated here because a referee reading the limitations alone must find it.
The composition baseline's 19 features span the 16-cell dinucleotide simplex, 15 degrees of
freedom, and those are exactly what the dinucleotide matcher controls: the GC arm controls 1 of
15, the dinucleotide arm 15 of 15. Had the model carried no information outside that simplex
both gains would collapse toward zero and the contrast would vanish, so the direction is
implied conditionally rather than guaranteed by algebra, which is what distinguishes it from
the composition-share quantity retracted below. But nobody should be persuaded by the sign. The
magnitude, its replication across cell lines, and the efficiency it buys are the claims.

**1e. Three architectures, not a survey.** This limitation used to read "one model class, and
the title should not outrun it", and R1g retired it: the contrast is measured for an
L2-penalised logistic on 4-mer counts, a 7,089-parameter CNN and a 19.7M-parameter fine-tuned
SpliceBERT, and it grows across all three. What remains is narrower and should be stated
plainly. Three architectures are not the space of sequence models. Nothing above 20M parameters
was tested, so the 100M-parameter LoRA-adapted models in the same family (RNA-FM, RNA-MSM) are
untested, and a model pretrained on multiple-sequence alignments could behave differently
because evolutionary context is not a property of the window's own composition. The ladder is
monotone over the range measured; whether it keeps rising, saturates or turns over beyond 20M
parameters **is not measured here**. The CNN step is also the weak rung, holding on 58 of 94
datasets, so the monotonicity claim rests mainly on the SpliceBERT end.

**1b. Negatives are not filtered for expression.** Region pools are built from all GENCODE v45
transcripts with no `gene_type` and no cell-line expression filter, and negatives are drawn
from the positive's chromosome rather than from within its transcript. A negative may therefore
sit in a gene not transcribed in K562 or HepG2. Combined with limitation 1, a substantial share
of the negative set is not RNA present in the cell, and the contrast is partly "RNA versus
not-RNA" rather than "bound versus unbound". ENCODE total RNA-seq exists for both lines; an
expression filter is the second revision we owe.

**1c. eCLIP peaks carry no significance threshold.** Every line of the ENCODE narrowPeak file
is ingested; the log10 p-value and log2 fold-change columns are never read, and no IDR or
fold-change filter is applied. The community standard since Van Nostrand 2020 is p ≤ 0.001 and
FC ≥ 4. Two spot-checked datasets already satisfy it (100% of peaks), so the practical effect
may be nil, but it is unfiltered by construction rather than by measurement. Note also that
`config/params.yaml` states the panel selects `preferred_default` experiments while
`src/rbp/data/encode.py` selects by replicate count, the code is right and the documentation
is wrong.

**2. Region-class assignment is asymmetric between arms, which matters for R1f.** Positives are classified by a
priority order over overlapping annotations; negatives are drawn from a merged pool of the same
region class. End to end, **6.1–8.7%** of negatives carry a region label the positive-side
classifier would not have given them. This is a second-order mismatch within an arm, not
between arms, and both arms share it, so it cannot generate the contrast. It does blur the
region split in R1f: a few percent of datasets may be assigned to the wrong dominant class,
which attenuates the CDS-versus-intron difference rather than creating it.

**3. Review process.** Every claim here has been through ten rounds of structured adversarial
critique with pre-registered acceptance criteria, and none of it has been through external peer
review. See `docs/59-the-council-and-the-correction.md`, which records the criteria, the
retractions, and a self-audit of whether those criteria were satisfied by construction.

---

## What the reader should take away

Three things any RBP variant-effect claim must clear before it means anything:

1. **A composition baseline under the negative-set protocol actually used.** Under GC matching,
   composition alone reaches 0.783 and the sequence model adds only **+0.0265**. Under
   dinucleotide matching the same model adds **+0.0662**. The protocol, not the model, decides
   most of what the number means.
2. **A positional baseline.** A leave-one-out prevalence rule over genomic blocks reaches
   **0.816** at 1 Mb and holds **0.755 to 0.817 within every phyloP decile**, so it is a second
   leak and not conservation in disguise.
3. **A conservation control read against a calibrated null, not against zero.** phyloP is
   strongly predictive here, and because logistic regression is not collapsible a *genuinely*
   independent covariate this strong would **amplify** the model's coefficient by roughly 65%.
   So "the coefficient barely moved when we adjusted for conservation" is evidence of
   substantial sharing, not of independence. Reporting an attenuation without simulating the
   null it should be compared against is the single easiest way to publish a backwards result,
   and this paper did exactly that before catching it (R4d).

**What this paper does not establish.** Whether the model beats conservation on ClinVar. The
AUROC ladder that appeared to settle it is circular, because ClinVar pathogenic assertions lean
on PP3 evidence from conservation-derived tools, and the coefficient analysis that replaced it
was computed wrongly and is retracted. R4 is unresolved and is labelled as such.
