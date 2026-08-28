# 60. The paper: every claim, every number, every limitation

**Manuscript source of truth.** Written 2026-08-27 after seven rounds of adversarial review.
If a number in the manuscript disagrees with this file, this file is wrong and must be fixed,
not the manuscript quietly edited.

**What the harness does and does not establish, stated precisely because the previous wording
here was false.** This header used to read "every number here is in a committed table and
gated by `scripts/verify.py`". It is not. `scripts/audit_manuscript.py` reports **45 numbers
in this document that trace to no committed table and no golden key** (see
`manuscript_orphans.csv`); most are aggregates over a subset that no whole-column aggregate
reproduces, but until each is traced the blanket claim is untrue. What is true:

- **150 numeric assertions** run in `scripts/verify.py` against `config/golden.yaml`, and the
  count of assertions that ran is itself asserted, so a gate cannot silently skip.
- Most of those are **regression gates**: they detect that a number changed, not that it was
  correctly derived, and several corruption attacks used to pass them. Call this a
  reproducibility harness, not verification.
- The exception is `scripts/recompute.py`, which **proves** rather than reproduces: 285
  published AUROCs rebuilt from committed per-example scores at max|difference| 2.2e-16, and
  it fails if the evidence is zeroed or deleted.

**Working title:** *What GC-matched negatives and ClinVar AUROC actually measure: trivial-baseline
calibration of RNA-binding-protein models across 94 ENCODE eCLIP datasets*

**Venue:** bioRxiv immediately; then *NAR Genomics & Bioinformatics* (Methods/benchmarking).

---

## The one-sentence claim

> Under negatives matched only on GC content, the protocol in general use, a 4-mer model's
> measured contribution over a mono+dinucleotide baseline is **+0.0265 AUROC**; matching on
> full dinucleotide composition raises it to **+0.0662** (difference **+0.0397** [+0.0336,
> +0.0458], larger in 88/94 datasets) while lowering the model's apparent AUROC by **0.1095**
> in **94/94** datasets. One fifth of that difference is arithmetic rather than biology,
> because a nested AUROC gain compresses against a higher baseline; net of compression the
> protocol effect is **+0.0313** [+0.0267, +0.0363], positive in 87/94. The standard protocol
> simultaneously inflates reported performance and conceals about half of the model's
> measurable contribution.

---

## R1: the protocol effect (primary result)

`cost_of_matching.csv`, n = 94 paired datasets. Figure **f1**.

| quantity | GC-matched | dinuc-matched |
|---|---|---|
| composition alone (19 features) | 0.7827 | 0.6280 |
| 4-mer score alone | 0.7981 | 0.6886 |
| composition + 4-mer score | 0.8092 | 0.6941 |
| **nested contribution of the score** | **+0.0265** [+0.0210, +0.0324] | **+0.0662** [+0.0559, +0.0765] |
| datasets where the score adds significantly | **80/94** | 82/94 |
| datasets where composition ≥ the k-mer model | **29/94** | n/a |

An earlier draft of this table put 0.7981 and 0.6886 on the "composition + score" row. Those
are the *standalone* 4-mer AUROCs. The nested values are 0.8092 and 0.6941, and the reported
contributions were always the nested ones, so the claim did not change; the row label was
wrong and is corrected here.

- cost of proper matching: **−0.1095** AUROC, **94/94 datasets fall**, paired Wilcoxon p < 1e-15
- contrast in nested contribution: **+0.0397** [+0.0336, +0.0458]
- **effect modification by dataset size is real: rho = +0.307, p = 0.0026** on the reported
  nested contrast. The paired design still rules out size *confounding* (both arms use the
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

| component of the +0.0397 contrast | value | share |
|---|---|---|
| attributable to AUROC compression | **+0.0083** [+0.0061, +0.0109] | 21% |
| **protocol effect, net of compression** | **+0.0313** [+0.0267, +0.0363] | **79%** |

The protocol effect is positive in **87/94** datasets and its interval excludes zero. The same
comparison computed directly on the unbounded d′ scale gives **+0.1290** [+0.1091, +0.1499],
larger in 87/94, which would be zero if compression were the whole story. Corrected for
compression, GC matching conceals **47.4%** of the model's measurable contribution, so the
honest phrase is "about half", not "two-thirds" (the uncorrected figure is 60%, and an earlier
draft said two-thirds, which was wrong on both counts).

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

So two of three scales agree with R1 and the third is measuring task difficulty rather than
incremental value. Both the reversal and its fingerprint are gated, so if the fingerprint ever
inverts the verifier fails and the reversal becomes real evidence against R1.

**Provenance note, stated because it matters.** Before this analysis the +0.0397 contrast
appeared in the manuscript, in no committed table, and under no golden key. It could not be
recomputed and it could not fail. It is now derived by `scripts/scale_check.py` from the two
rehearsal arms and asserted by `r1_scale_check`. Its interval, recomputed with a paired
dataset-level bootstrap (2000 draws, seed 0), is [+0.0336, +0.0458]; the manuscript previously
printed [+0.0334, +0.0461] from an unrecorded computation.

## R1c: the strand artifact, measured against a matched placebo

`strand_placebo.csv`, `scripts/strand_placebo.py`, 16 gated checks, n = 40 datasets.

**The wound.** `annotation.py:126` drops region strand by design ("a window's strand comes from
its peak"), which is right for positives and wrong for negatives, so `negatives.py:328` gives
each negative the *positive's* strand. Only **55.6%** of negatives sit on the strand their own
gene is transcribed from; positives are 100% true sense. Direction is therefore a cue that
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
| placebo, same n dropped at random | +0.0346 [+0.0258, +0.0444] |
| placebo, matched on the retained region marginals | +0.0322 [+0.0234, +0.0418] |

Decomposing the −0.0091 that restriction alone reports: **−0.0032** is the cost of discarding
pairs at all, **−0.0024** [−0.0048, −0.0003] is locus mix, and the remainder,
**−0.0036** [−0.0071, +0.0001], is strand. The strand-corrected contrast is
**+0.0342** [+0.0253, +0.0442] and **90.6%** of the effect survives.

**Correction, and it goes the paper's way.** An earlier version of this section used the
unstratified placebo, reported the excess as −0.0059 with an interval clear of zero, and said
"the artifact is real". Once the placebo is matched on region the interval is [−0.0071,
+0.0001] and that statement is withdrawn: the strand artifact is **small and not
distinguishable from zero**, and only a bound is claimed. The component that *is*
distinguishable from zero is the locus mix, which is why the plain placebo was not a valid
counterfactual in the first place.

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

`r1_robustness.csv`, `scripts/r1_robustness.py`, 11 gated checks.

R1 concedes that the *sign* of the contrast is implied by the design, so only the magnitude is
informative. This section is what the magnitude is worth. Neither result below is implied by
the design.

**It replicates out of sample.** Fifteen proteins were assayed in both HepG2 and K562. Those
are separate eCLIP experiments with separately drawn negatives, so the contrast measured in one
line is an out-of-sample prediction of the other. Pearson **r = +0.909** [+0.812, +0.972],
p = 2.6e-06; Spearman +0.932. Mean absolute difference between lines **0.0151**, against a
between-dataset standard deviation of 0.0318.

Crucially the contrast replicates **better than either arm it is built from** (GC-arm gain
alone r = +0.518, dinucleotide-arm gain alone r = +0.813). The difference is a more stable
property of the protein than either measurement, which is the opposite of what noise would do.
The design guarantees the sign in each cell line independently and guarantees nothing about
whether the magnitudes agree, so this is information about magnitude alone.

**It buys statistical efficiency, which is what a benchmark builder acts on.** The nested gain
rises 2.50x under dinucleotide matching, but its standard error rises only 2.18x, so the
per-dataset signal-to-noise z = gain / SE rises **1.31x** (mean 9.97 to 13.11, median 7.01 to
10.13), higher in **83/94**, paired Wilcoxon p = 3.1e-14. Because z grows as the square root of
sample size, the same contribution reaches the same confidence with **58%** of the labelled
windows. And composition alone beats the sequence model in **29/94** datasets under GC matching
against **14/94** under dinucleotide matching: the harder protocol makes the model look better
on that comparison, not worse.

**This is the practical recommendation, and it does not depend on either arm being "the
truth".** A benchmark built on dinucleotide-matched negatives reports a lower headline AUROC
and a more precisely measured, more reproducible, and roughly twice as large estimate of what
the model actually contributes.

## R1e: the contrast rebuilt from sequence, and independent of k

`k_sweep.csv`, `scripts/k_sweep.py`, 11 gated checks, all 94 datasets, both arms.

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
| **R6**, the substitution-spectrum baseline | Its interval, [−0.0113, +0.0284], is wider than half the effect. |
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

**1. Negative strand assignment, and the primary result has no strand control.** Sampled
negatives inherit the *positive's* strand while their location is chosen to match composition,
so strand is effectively a coin flip. Audited on 40 of 95 datasets: **55.2%** of negatives carry
the strand of the gene they sit in, but that is the fraction among *unambiguous* windows only.
A further **14.0%** lie in genes transcribed on both strands and are excluded from the
denominator, so the demonstrably-sense fraction is **47.4%**. 0.0% are intergenic. Roughly half
the negative set is therefore antisense sequence no transcript produces, while every positive
is true sense RNA. This inflates all absolute AUROCs in R1.

**An earlier draft defended this with the wrong statistic and the defence is withdrawn.** The
numbers offered (contrast +0.2643 → +0.2787 on sense-only negatives; +0.2485 vs +0.2644 across
panel halves) are the *composition-share contrast* (the quantity retracted below as an
algebraic identity) and one of them appeared in no committed table. R1's reported contrast is
the **+0.0397 nested gain**, and no strand-restricted recomputation of it exists. What survives
from that analysis is narrower but real: the directional cue is detectable and
SpliceBERT-specific (it ranks sense above antisense negatives at 0.576 against the CNN's 0.521,
paired p=3.2e-6), so a pretrained RNA model can read strand. Whether that changes the +0.0397
is **untested**. Regenerating the negatives with strand-correct sampling requires only the
composition and k-mer arms, which are CPU-only, and is the first revision we owe.

*A previous draft also asserted that this convention is shared with the matched-negative
samplers this work benchmarks against. That sentence is deleted: we have no citation for it,
and it is probably false: RNA samplers such as GraphProt draw unbound windows from annotated
transcripts on the transcribed strand, which is transcript space, not genome space.*

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

**1d. Conservation is partly a re-read of the labels.** ClinVar pathogenic assertions for
noncoding SNVs rest substantially on PVS1 and PP3, and PP3 is operationalised through
CADD/REVEL/GERP, all conservation-derived; benign assertions lean on BA1/BS1/BP4, i.e. allele
frequency. In this study's own data benign variants have mean phyloP **0.022** (median −0.046,
1.1% above phyloP 7) against pathogenic mean **5.06** (38.4% above phyloP 7). phyloP is
therefore not an independent predictor of these labels but a partial proxy for the evidence
that generated them, and its 0.892 must be read as a **ceiling on this benchmark**, never as a
biological statement. The ClinVar sections that rested on it are cut, and this
paragraph is retained only to record why: a benchmark whose labels are partly derived from the
competitor being benchmarked cannot rank the two, however carefully the ranking is computed.

**2. Pretraining exposure.** SpliceBERT is pretrained on human pre-mRNA. Every window in this
study (positive and negative) lies inside an annotated gene (0.0% intergenic), so held-out
chromosomes are plausibly represented in its pretraining corpus. It is therefore the only arm
carrying information from the test folds, and its AUROC is an upper bound relative to the three
arms trained from scratch. This cannot be quantified without the pretraining corpus, applies
identically to every published genomic language model, and is one reason the primary result
(R1) is stated without reference to SpliceBERT.

**3. ClinVar assertion quality.** No review-status filter is applied: `CLNREVSTAT` is not
parsed, so zero-star and single-submitter assertions are included. Exact CLNSIG matching does
exclude Uncertain-significance and Conflicting records. A ≥1-star restriction would reduce the
variant set and is the obvious first robustness request.

**4. Region-class assignment is asymmetric between arms.** Positives are classified by a
priority order over overlapping annotations; negatives are drawn from a merged pool of the same
region class. End to end, **6.1–8.7%** of negatives carry a region label the positive-side
classifier would not have given them. This is a second-order mismatch within an arm, not
between arms, and both arms share it.

**5. Cross-dataset variant sharing.** 68% of a target's variants appear in at least one other
powered dataset. Inference clusters on 1-Mb genomic blocks, and a cluster bootstrap over
sharing components gives [0.0215, 0.0998] against [0.0269, 0.0993] iid; the intervals overlap
substantially, so the sharing is real but not materially inflating.

**6. The specificity estimate is not corrected for errors in a covariate.** Donor quality is a
measured binding AUROC with its own sampling error, so regression dilution biases the R5
intercept toward the unadjusted mean. A model-free local estimate on near-neutral donor pairs
(|Δquality| < 0.05, |Δlog size| < 0.25) gives **+0.0774** [+0.0360, +0.1163], 14/18, p=0.0034,
agreeing with the regression intercept to three decimals, so the bias appears small, but it is
not formally bounded.

**7. Review process.** Every claim here has been through six rounds of structured adversarial
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

And one design requirement: **a wrong-protein control must match donors on model capacity, and
cannot detect specificity at all below ~20 positive variants** (with the gap still mid-slope
there, and its plateau only at ≥45).

**What this paper does not establish.** Whether the model beats conservation on ClinVar. The
AUROC ladder that appeared to settle it is circular, because ClinVar pathogenic assertions lean
on PP3 evidence from conservation-derived tools, and the coefficient analysis that replaced it
was computed wrongly and is retracted. R4 is unresolved and is labelled as such.
