# 59. The council, the correction, and the experiment that ended it

**Covers 2026-08-26 (evening) to 2026-08-27.** Doc 58 ends with the run complete, twenty-four
bugs catalogued, and R4 "complete, and it took three corrections to become honest". Doc 55
ends at Bug 31. This file picks up from there and covers what happened when the work was
handed to four rounds of adversarial review, what those rounds destroyed, the $4 experiment
that settled the central question, and the state the project is now in.

**Read this if:** you want to know why the paper's headline changed, why R4 is no longer a
positive claim, and why the verifier's "80/80 passed" is both better and weaker than it looks.

**One-line summary.** Three independent reviewers found that the wrong-protein control was
confounded by donor model quality; a re-run with five quality-spanning donors per target
removed the confound and showed the effect is real but only in well-powered datasets; the
paper's headline moved to a finding nobody had noticed we already had; and the single most
striking claim turned out to be published in 2015.

---

# Act 0. The state this file starts from

At the end of doc 58 the project stood at:

| thing | value |
|---|---|
| golden checks | 71/71 passing, fetched from GCS |
| unit tests | 576 passing, 6 skipped |
| result tables | 22 |
| figures | 6 (png + pdf) |
| bugs catalogued | 31 |
| money spent | ~$30 Modal free credit + GCP credits |
| manuscript | not started |

And the claim set was:

- **R1** — composition reproduces 94.8% of the AUROC gain under GC-matched negatives, 67.8%
  under dinucleotide-matched; difference 27.0 points.
- **R2** — composition 0.628 < k-mer 0.688 < CNN 0.706 < SpliceBERT 0.809. *Demoted to a
  methods table* (Horlacher 2023 got there first).
- **R3** — ISM locality. *Cut* (architecture-confounded by construction).
- **R4** — the wrong-protein control. Paired gap **+0.0645** on 44 powered datasets, 33/44,
  p=3.9e-04. Clean-donor stratum **+0.121**, 16/17, p=1.5e-04. **This was the paper's one
  novel positive claim.**
- **R5** — a trivial 1-Mb positional-prevalence rule (0.818) and phyloP (0.892) both beat the
  model (0.755). **This was believed to be the striking result.**

Every one of those five lines is different now. Here is how.

---

# Act I. Why a fourth round, and how it was made to terminate

## The problem with the first three rounds

Rounds 1-3 each found real problems, and each time the fix produced a new artifact that the
next round could attack. From the outside this is indistinguishable from an infinite loop, and
that is exactly how it felt. It was not actually a loop — the rounds were walking down a
ladder of abstraction:

| round | attacked | found |
|---|---|---|
| 1 | the claims | R4's pooled AUROC was Simpson-inflated (0.829 → 0.755) |
| 2 | the framing | R2 and R3 were already published or confounded; cut |
| 3 | the code | no CI, a dead golden key, an `abs()` that shrank the effect |

Round 3 found **zero wrong numbers**. That is the signature of convergence, not a treadmill.
But there was no *stopping rule*, and an adversarial reviewer asked "find problems" will always
find something. The process could not terminate on its own.

## The termination mechanism

The fix was to change what the council was asked for. Round 4 required every advisor to
**pre-commit, in writing, before seeing any new data**, to:

- **B1** — the exact numerical criterion on which they would accept the effect exists
- **B2** — the exact criterion on which they would accept a null
- **B3** — what would be ambiguous, and what to *write* in that case (not "run more")
- **B4** — any confounder not yet on the list

Then the experiment ran, and each advisor applied their own stated criterion mechanically.
This works because the goalposts are fixed before the ball is kicked. There is nothing left to
argue about afterwards, and an advisor who moves their own threshold is visibly doing so.

All four advisors who pre-committed landed on the same branch. That is the reason this file
exists and a doc 60 does not.

## The remits, chosen not to overlap

Rounds 1-3 had attacked the claims, the framing, the statistics and the code. Round 4 was
pointed at what nobody had touched:

| advisor | remit | why it was new |
|---|---|---|
| A | **data provenance** — ENCODE selection, negative construction, fold leakage, ClinVar processing | nobody had questioned the *data*, only the analysis of it |
| B | **independent reproduction** — recompute every headline from scratch, then try to break the verifier | nobody had recomputed a single number without using the pipeline's own code |
| C | **prior art** — hunt for papers that already published each surviving claim | round 2's literature pass was shallow |
| D | **whole-artifact coherence** — does the MLOps story and the paper form one thing | nobody had asked whether the pieces fit |

A fifth agent independently audited whether the *documentation* claims matched the code.

---

# Act II. The collapse of R4

This is the longest section because it is the most important thing that happened, and because
the sequence of "the claim is dead / no it isn't / yes it is / actually it's conditional" is
the part that is easy to get wrong when retelling.

## II.1 The finding: the donor was also the worse model

Three advisors, independently and by different routes, converged on the same defect.

The wrong-protein control scores each dataset's ClinVar variants with a *different* protein's
fine-tuned head. The donor was chosen in `scripts/variant_splicebert.py` as:

```python
w = man.iloc[(index + mismatch) % len(man)]      # mismatch = 47
```

One donor per target, one deterministic offset. The measured consequences, on the 44 powered
datasets:

| statistic | value |
|---|---|
| donor trained on fewer pairs than target | 28/44, Wilcoxon p=0.023 |
| donor's own binding AUROC vs target's | **0.802 vs 0.850**, p=0.018 |
| gap vs log10(donor pairs) | **rho = −0.533, p=1.9e-04** |
| ↳ partialling out power | −0.497, p=6.0e-04 |
| ↳ partialling out power **and** target size | −0.493, p=6.7e-04 |
| the same for *target* size, partialling power | −0.017, p=0.91 (nothing) |
| the floor itself vs donor pairs (partial on power) | **+0.481, p=9.4e-04** |

Read the last two rows together. It is the **donor's** training volume that predicts the gap,
not the target's. And the wrong-protein floor — the quantity the whole control assumes is
protein-agnostic — tracks how much data the donor saw.

So "the right protein's head beats the wrong protein's head" partly meant **"my head was
trained on more data than that head."**

## II.2 The placebo, which is what actually killed it

The published claim was stratified on donor-target *co-binding overlap*: the 17 donors with
negligible overlap gave +0.121, 16/17. An advisor split the same 44 datasets on **donor size**
instead — a variable with nothing to do with co-binding:

| split of the same 44 datasets | n | gap | wins | p |
|---|---|---|---|---|
| **published** — clean donor by peak overlap | 17 | +0.1210 | 16/17 | 2.7e-04 |
| **placebo** — the 17 *smallest* donors | 17 | **+0.1362** | 16/17 | 2.7e-04 |

The placebo reproduces the flagship result **better** than the real stratifier. That is what a
confounded stratification looks like from the outside.

And the decisive cut:

| stratum | n | gap | wins | p |
|---|---|---|---|---|
| donor's own binding model ≥ target's | 18 | **−0.0025** | 9/18 | **1.00** |
| donor's model weaker than target's | 26 | +0.1109 | 24/26 | 1e-05 |

Where the donor was not the weaker model, there was nothing.

## II.3 Why the contamination metric was measuring the wrong thing

The contamination test that had been added to *defend* R4 was itself broken:

```python
"shared_frac": len(t & d) / len(t),     # normalised by the TARGET only
```

Because the denominator is the target's variant set, a donor with a small variant set *cannot*
score high. Measured:

- `shared_frac` vs donor variant-set size: **rho = +0.696, p=1.6e-07**
- "clean" donors: median **40** assigned variants. "dirty" donors: median **320**
- the flagship claim (floor vs overlap, rho=+0.299, p=0.049), **controlling for donor set
  size: rho = −0.240, p=0.12** — reverses and dies

The size-symmetric measure, `jaccard`, was already computed in the same table and was never
used for the stratification. Using it, the direction survives but the correlation claim does
not.

## II.4 The rescue attempt, and why it failed

Stratifying throws away data and re-confounds on power. So rather than stratify, the gap was
regressed on the confounders directly and the **intercept** read — the estimated gap when
donor and target are equally good and equally sized. 5000-draw percentile bootstrap:

| adjusted for | intercept | 95% CI | p |
|---|---|---|---|
| donor quality | +0.0442 | [+0.0006, +0.0845] | 0.048 |
| donor size | +0.0399 | [−0.0002, +0.0747] | 0.051 |
| **both** | **+0.0372** | **[−0.0055, +0.0764]** | **0.087** |

Roughly half the published +0.0645, and not significant with both confounders in. But still
positive in every specification — so "dead" was too strong.

**Then the rescue was itself killed.** The statistician ran the identical regression on all 82
usable datasets instead of the 44 powered ones:

| adjusted for | powered 44 | all usable 82 |
|---|---|---|
| quality + size | +0.0372 [−0.0055, +0.0764] | **−0.0142 [−0.0515, +0.0214]**, p=0.46 |

**The sign of the adjusted intercept is decided by the stratum boundary.** You cannot appeal
to a regression to escape a stratification when the regression's answer flips with the
stratification. Reproduced independently and confirmed exactly.

That is what "not identified" means. R4 as published was not salvageable from existing data.

## II.5 One advisor overstated, and it mattered

The genomics reviewer reported the confound-free half as a flat null and concluded the
question "cannot currently be answered." That was too strong, and the correction is worth
recording because it shaped the experiment:

The donor≥target stratum is confounded **in the opposite direction**. Its median power is
**54.5** pathogenic variants against **151.5** in the complement (MWU p=2.2e-04), and the gap
is known to grow with power. At a higher power threshold the same stratum gives +0.059, 6/8.

So the "clean" stratum traded one confound for another. Neither stratification was
trustworthy. That is precisely why the fix had to be a *new draw*, not a re-analysis.

## II.6 A factual dispute between advisors, resolved

Two advisors disagreed about *why* donors were systematically weaker.

- One said the manifest is sorted by pair rank, so `+47` deliberately picks a differently
  sized dataset. The code comment said the same.
- One said the manifest is alphabetical, and measured target/donor pair-rank correlation at
  −0.079.

Settled by reconstructing both orderings and applying `+47`: **alphabetical predicts the
observed donor 82/82; pair rank predicts 1/82.** The manifest is built by
`a.groupby(["protein","cell"])`, and pandas sorts groupby keys.

Two consequences:

1. Two code comments were **false statements**, not merely stale — `variant_splicebert.py`
   ("pair-rank-sorted manifest") and `modal_variants.py` ("sorted by pair rank, so a dataset
   is paired with one of a very different size"). Both now deleted rather than corrected.
2. **Had the comment been true, the confound would have been worse.** Under a genuinely
   pair-rank-sorted manifest with `+47`, target/donor size correlation would be −0.497
   (p=2.1e-06), hard-wiring the anti-correlation in donor capacity. Observed is −0.079. The
   accident of alphabetical ordering made the confound random rather than systematic.

---

# Act III. The experiment

## III.1 Design decisions, and the one that was rejected

The obvious fix is to **match** donors to targets on size and quality. It was tried on paper
and rejected for two independent reasons.

**Reason 1 — it destroys the estimand.** Matching removes the variance in donor quality, and
donor quality is precisely the regressor whose slope identifies the intercept. Match on it and
there is nothing left to adjust for.

**Reason 2 — the pool is empty exactly where it matters.** Measured before spending anything
(`scratchpad/pool.py`), over the 95-dataset pool:

| constraint | median eligible donors | targets with ≥5 | targets with 0 |
|---|---|---|---|
| \|Δlog10 pairs\| ≤ 0.15 **and** \|ΔAUROC\| ≤ 0.05 | 8.0 | **36/44** | **3** |
| \|Δlog10 pairs\| ≤ 0.15 only | 17.5 | 41/44 | 0 |
| \|Δlog10 pairs\| ≤ 0.25 only | 29.0 | 43/44 | 0 |
| \|Δlog10 pairs\| ≤ 0.40 only | 48.0 | 44/44 | 0 |

Targets are systematically the *strong* datasets, so a tight match starves the pool for the
datasets carrying the claim. This check cost $0 and ran before any GPU time — it is the
"measure before you commit" discipline applied to an experimental design rather than a bill.

**The design chosen: span, don't match.** Five donors per target, drawn to cover the donor
quality range, with a hard constraint of ≥2 donors *stronger* than the target where the pool
allows. Spanning gives support on both sides of zero, so the estimand is an **intercept**
rather than a subset mean.

**Screening changed too:** donors screened on `jaccard` ≤ 0.02, not the broken `shared_frac`.

## III.2 The realised draw

`scripts/donor_draw.py` (168 lines). Output `results/tables/donor_tasks.tsv`, uploaded to
`gs://rbp-repro-2026-derived/variants/donor_tasks.tsv`.

```
475 tasks, 95 targets, 5-5 donors each
donor advantage (donor - target binding AUROC): min -0.296 median +0.036 max +0.296
  donors STRONGER than target: 328/475 (94 targets have >=1)
log10 donor/target pairs: min -1.69 median -0.06 max +1.76
jaccard: median 0.0000 max 0.0199
targets with >=2 stronger donors: 92/95
```

Compare with the old arm, where donors were systematically *weaker*. **328 of 475 donors are
now stronger than their target.** That inversion is the whole point.

## III.3 The code

Three changes, deliberately small, because the arms must stay comparable.

**`scripts/variant_splicebert.py`** — a `--donor-task` path that reads the donor manifest and
names target and donor explicitly instead of deriving the donor from an offset:

```python
if donor_task >= 0:
    dm = pd.read_csv(io.StringIO(bucket.blob(DONORS).download_as_text()), sep="\t")
    d = dm.iloc[donor_task]
    index = int(d.target_idx)
    r, w = man.iloc[index], man.iloc[int(d.donor_idx)]
    out = f"variants/scores_md/{cell}_{prot}__{wcell}_{wprot}.csv"
```

The scoring body was factored out into `_score(...)` so that all three arms — matched,
single-donor mismatched, multi-donor — run *literally the same code*. Duplicating a scoring
body per arm is how arms silently drift apart.

**`cloud/modal/modal_variants.py`** — a `multidonor_sweep` entrypoint. One non-obvious line:

```python
rcs = list(task.starmap([(0, force, 0, i) for i in range(n)]))
```

`starmap`, not `map`. `map` varies the *first* positional argument (`idx`), and the thing that
has to vary here is `donor_task`. This was caught before launch by reasoning about the
signature; `map` would have run task 0 with donor 0, 475 times, and written one file.

**`MAX_CONTAINERS` 10 → 20**, with the arithmetic written next to it: this is the T4 inference
arm, so the ceiling is 20 × $0.59 = $11.80/h, and the 30-minute per-task timeout bounds a
total hang at 20 × 0.5 × $0.59 = $5.90. The A10G *training* sweep keeps its own lower cap in
`modal_sweep.py`, where the same concurrency would cost $22/h. Raising this one halves wall
time and changes the bill by the price of ten extra cold starts, because Modal bills
container-seconds, not containers.

## III.4 The run

A **one-task probe first**, per this project's own rule:

```
[04:06:05] AATF:K562 weights=MTPAP:K562 54 variants on cuda in 3s
1/1 ok
```

Then the full sweep, detached:

| | |
|---|---|
| tasks | **474/475 ok** |
| scored rows | 164,176 |
| targets | 95 |
| distinct donors used | 76 |
| usable (target, donor) AUROC pairs | 409 |
| wall time | ~20 min |
| **cost** | **~$4** |

`scripts/scores_md/` is under `variants/`, which `rbp-modal` is permitted to write — checked
against the IAM table before launching, not after.

## III.5 The analysis

`scripts/multidonor_analysis.py` (230 lines). Three estimators, because one would have been a
choice and three are a robustness check:

1. **Target-clustered bootstrap intercept.** Resamples *targets*, not rows — five rows share a
   target and are not independent.
2. **Within-target adjusted.** Centres gap and donor advantage within each target, so target
   identity cannot contribute at all. Every between-target difference cancels.
3. **Unadjusted mean gap**, targets as the unit.

Both panels are always reported. The old arm's intercept was +0.037 on the powered 44 and
−0.014 on all 82 — reporting one panel and not the other is how that stayed hidden.

## III.6 The result

**Powered stratum (n_pathogenic ≥ 20): 44 targets, 219 pairs**

| estimator | value | 95% CI | p |
|---|---|---|---|
| mean gap, unadjusted | **+0.0880** | — | 4.6e-08 (39/44) |
| intercept given advantage | +0.0879 | [+0.0639, +0.1106] | ~0 |
| intercept given advantage+size | **+0.0784** | [+0.0514, +0.1038] | ~0 |
| intercept given advantage+size+**power** | +0.0236 | [−0.0120, +0.0593] | 0.198 |
| within-target adjusted | +0.0881 | [+0.0630, +0.1112] | ~0 |
| donors stronger than target only | **+0.0920** | — | 1.3e-07 (39/44) |

**All usable: 82 targets, 409 pairs**

| estimator | value | 95% CI | p |
|---|---|---|---|
| mean gap, unadjusted | +0.0124 | — | 0.104 (56/82) |
| intercept given advantage+size | +0.0068 | [−0.0278, +0.0399] | 0.69 |
| donors stronger than target only | +0.0090 | — | 0.198 (53/81) |

**The confound is verifiably gone:**

| diagnostic | old | new |
|---|---|---|
| gap vs donor advantage | rho **−0.533**, p=1.9e-04 | rho **−0.028**, p=0.57 |
| donors stronger than target | −0.0025, 9/18, **p=1.00** | **+0.0920, 39/44, p=1.3e-07** |

The old design's decisive null was a one-donor artifact. With 141 stronger-donor pairs across
44 targets instead of 18 single draws, that stratum became the *strongest* evidence.

## III.7 Why the two panels disagree — mechanism, measured

The low-power stratum is still significantly negative: n=38, **−0.0750, 17/38, p=0.025**. So
it survived the fix and is not a donor artifact. The mechanism was measured rather than
assumed:

| stratum | matched arm | wrong-protein floor |
|---|---|---|
| powered 44 (median 90 pathogenic) | 0.7556 | 0.6675 |
| low-power 38 (median 9 pathogenic) | **0.5593** | 0.6343 |
| change | **−0.1963** | −0.0332 |

AUROC vs log10(power), across all 409 pairs:

- matched arm: **rho = +0.631, p=7.6e-47**
- wrong-protein floor: **rho = +0.091, p=0.066** — flat, as a generic signal should be

Gap by power bin: 1-5 pathogenic **−0.107** (sd 0.352) | 5-20 **−0.065** | 20-100 **+0.056** |
100+ **+0.133** (sd 0.074).

**It is the matched arm that collapses, not the floor that rises.** Below ~20 pathogenic
variants the target's own head sits at 0.559 — near chance — because low-power datasets are
also small datasets whose heads are undertrained (target binding AUROC 0.791 vs 0.838,
p=0.0063). A near-constant |δ| yields AUROC → 0.5, while the donor head, often trained on more
data, is still a functional scorer. That is a **degenerate-scorer comparison**, not negative
specificity.

Also relevant: 10 pairs across 2 targets have <20 *total* variants (mean gap −0.236).
Excluding them, the all-82 result becomes +0.0187, 56/80, p=0.040.

## III.8 The verdict

Every advisor applied their own pre-committed criterion. All of them fired **B3 — ambiguous /
conditional**, and all for the same reason: the intercept excludes zero on the powered 44 and
includes it on all 82.

| advisor | their B1 threshold | outcome |
|---|---|---|
| B (reproduction) | intercept ≥ +0.02, CI excluding 0 on **both** panels | fails on the 82 |
| C (prior art) | CI lower bound > +0.02 on **both** panels | fails on the 82 |
| D (coherence) | intercept ≥ +0.030 on 44 **and** CI excluding 0 on 82 | fails on the 82 |

Nobody moved their threshold after seeing the data. Two advisors explicitly noted they were
not softening. One noted their own B4 concerns (no donor random effect across 76 donors used
~5 times each; no attenuation correction on a regressor measured with error) remain
unaddressed and inflate confidence in the +0.0784.

**So R4 is a Results subsection reporting both panels, and it is not in the abstract.**

The genuinely novel piece is not "specificity exists" — the control logic is owned by Adebayo
et al. 2018 (*Sanity checks for saliency maps*) and Hooker et al. 2019 (ROAR). It is the
**detection threshold**: the generic floor is power-invariant while the specific signal is not,
so specificity claims below roughly 20 positives are unmeasurable by construction. No prior
art was found for that calibration.

---

# Act IV. The strand bug — the most dangerous finding, and the test that saved the headline

## IV.1 Discovery

The data-provenance advisor read the negative-sampling code and found:

```python
# src/rbp/data/annotation.py:127
# "Strand is deliberately dropped. A window's strand comes from its peak, so the
#  region's own strand is never needed for classification or negative matching."

# src/rbp/data/negatives.py:328-330
"strand": p["strand"],                          # the POSITIVE's strand
"seq_rna": to_rna(neg["seq_dna"], p["strand"]), # revcomp decided by the POSITIVE
```

A negative's genomic location is chosen to match *composition*, not strand. Whether it lands
on a gene transcribed in the assigned direction is a coin flip.

Measured independently across 40 datasets (`scripts/strand_audit.py`, built for this):

```
negatives whose ASSIGNED strand matches their own gene: mean 55.2% (range 43.3%-61.5%, sd 0.038)
negatives in no annotated gene:        0.0%
negatives in genes on BOTH strands:    14.0%
```

So **~45% of negatives are antisense** — sequence no transcript in the cell produces — while
100% of positives are true sense RNA. The negatives are all inside genes (0% intergenic), so
"antisense" is a meaningful category here, not an artifact of the measurement.

## IV.2 Why this was the most dangerous thing found

Sense/antisense is a **non-compositional** cue. Dinucleotide matching is performed on forward
DNA and `revcomp` is applied to *both* members of a pair, so the composition match survives
into RNA space — but the directional cue does not. Real transcripts have directional features
(splice motifs, polypyrimidine tracts, polyA signals). Read backwards, those vanish.

A model pretrained on pre-mRNA can read direction. A mononucleotide+dinucleotide model largely
cannot.

**That is a complete alternative explanation for the new headline.** "Composition reproduces
only 41% of SpliceBERT's skill" could mean *SpliceBERT is not finding more motif, it is
noticing which negatives are backwards.*

## IV.3 The test, and why it was free

The GTF used to build the annotation (`gencode.v45`) was already on disk, and the processed
datasets with negative coordinates were in GCS. So the test needed no GPU and no retraining:

1. Build a gene-strand interval index from the GTF (63,187 genes, 25 contigs).
2. For each negative, find overlapping genes and compare their strand to the assigned strand.
3. Per-dataset **sense fraction**.
4. Correlate with the per-dataset shares.

**The key insight that makes this decisive:** the claim is a **contrast** between two shares,
not a level. For the artifact to *create* a contrast, it must act **differentially** on the
two models. So the test is not "does the artifact matter" (it does) but "does it matter more
for SpliceBERT than for the k-mer model."

## IV.4 The result: the bug is real, the headline survives

```
frac_sense vs THE CLAIM: k-mer share minus SpliceBERT share   rho +0.2253  p=0.162
frac_sense vs k-mer share                                     rho +0.4053  p=0.009
frac_sense vs SpliceBERT share                                rho +0.5167  p=0.001
frac_sense vs composition AUROC                               rho +0.4173  p=0.007
```

The artifact moves **both** shares hard and in the same direction. It does **not** move the
contrast (p=0.16, not significant).

And the model-free version, splitting the panel at the median sense fraction:

| half | k-mer share | SpliceBERT share | **contrast** |
|---|---|---|---|
| antisense-rich (n=20) | 0.5779 | 0.3294 | **+0.2485** |
| antisense-poor (n=20) | 0.7846 | 0.5201 | **+0.2644** |

The objection predicts a *much larger* contrast where the artifact is strongest. The data give
a slightly **smaller** one. The objection loses.

## IV.5 What this does not settle, stated plainly

- The test used **40 of 95** datasets.
- The sense fraction varies narrowly (sd 0.038, range 0.43-0.62), so a correlation across
  datasets has **limited power** to detect a differential effect. The stratified contrast is
  the more trustworthy half of the evidence, because it is a direct comparison rather than a
  regression on a barely-varying regressor.
- The bug **does** inflate every absolute number in R1. It is a real limitation and belongs in
  the paper as one.
- Fixing the sampler and re-running is the real answer. This bounds the damage in the meantime.

All of it is now gated in `config/golden.yaml` under `strand_audit`, including
`contrast_rho_p_min: 0.05` — an assertion that the correlation must stay **non-significant**.
If the artifact ever starts predicting the contrast, the build fails.

> **CORRECTION, 2026-08-27.** The paragraph above was false when written. The `strand_audit`
> block existed in `golden.yaml` and `grep strand scripts/verify.py` returned **nothing** — the
> nine keys were written down, not gated. A reviewer caught it, and it is the *third*
> occurrence of this exact bug: `integrity.min_tests_passing` sat unread while the suite grew
> from 480 to 576; `r1_headline_is_gc_share_only` was added to forbid a headline and, unread,
> failed to stop that headline being promoted two weeks later; and then the commit that wired
> up those 26 keys introduced 9 more that nothing read, plus this sentence asserting they were
> covered.
>
> `verify_strand_audit` now exists and all nine keys are read. More to the point, the defence
> is no longer a person remembering to check: `tests/unit/test_golden_keys_are_read.py` parses
> `golden.yaml`, walks every leaf key, and **fails the build** if any key is referenced
> nowhere. Exemptions are a named dict with a written reason per entry, so adding one is a
> visible act in a diff.
>
> The lesson is not "be more careful." It is that a claim about coverage is itself a claim, and
> claims need gates. This document asserting a gate exists is exactly as trustworthy as a
> config asserting a number is checked — which is to say, not at all, until something executes.

---

# Act V. The novelty demolition

The prior-art advisor was asked to search hard, and delivered the single most consequential
reframing.

## V.1 Verdicts

| claim | verdict | prior art |
|---|---|---|
| composition share is **model-dependent** | **INCREMENTAL but the one new number** | negative-set matching itself is old (Lee/Karchin/Beer 2011 kmer-SVM; Ghandi/Lee/Beer 2014 gkm-SVM; dinucleotide-preserving shuffles, Bailey 2011). Whalen/Schreiber/Noble/Pollard 2022 *Nat Rev Genet* names GC confounding as a canonical pitfall. **No prior art found for the share expressed as a fraction of gain over chance and shown to differ by model class.** |
| trivial positional rule and phyloP beat the model | **ALREADY PUBLISHED** | Grimm et al. 2015 *Human Mutation* — "type 2 circularity": a predictor knowing only *gene identity* scores high, and benchmark rankings invert under gene-disjoint evaluation. Schreiber/Singh/Bilmes/Noble 2020 *Genome Biology* — the *average activity* baseline at a genomic position beats deep models that appear to learn cell-type specificity. Plus Livesey & Marsh, ProteinGym (Notin 2023), BEND (Marin 2024), Tang & Koo, Sasse/Mostafavi 2023, Karollus & Gagneur 2023. |
| wrong-protein control must be matched on capacity | **ALREADY PUBLISHED as principle** | Adebayo et al. 2018 NeurIPS (model/label randomization sanity checks); Hooker et al. 2019 (ROAR); Raghu et al. 2019 (Transfusion). |

## V.2 What this changed

**The 1-Mb baseline result is a restatement, in a new tissue, of a conclusion held since
2015.** It was believed to be the striking finding. It is a *supporting null* — and running it
before a reviewer did is still a genuine credit, just not a contribution.

**The centre of gravity moved to R1's model-dependence**, which was sitting in the tables
unnoticed. Note its direction: it is a *pro*-language-model result, which is rarer than the
anti-LM literature and therefore more interesting.

## V.3 The number that became the headline

R1's share had always been computed against the **k-mer** model, because that is what the
rehearsal arm trains. It had been described in summaries as a fraction of "the four-model
gain," which is wrong. The repo was always honest — `figures.py:158` sets the f1 y-axis to
literally "k-mer model AUROC" — the error was in the prose.

Computed properly, same estimator (ratio of means, dataset bootstrap, 2000 draws), n=95:

| model | AUROC | composition share | 95% CI |
|---|---|---|---|
| k-mer | 0.6875 | 0.6818 | [0.6249, 0.7362] |
| CNN | 0.7063 | 0.6197 | [0.5649, 0.6748] |
| **SpliceBERT** | 0.8091 | **0.4137** | [0.3679, 0.4580] |

**Contrast: +0.2681 [+0.1970, +0.3383]**, excludes zero.

Two things make this trustworthy:

1. **Cross-check.** The k-mer row is the same quantity as the published dinuc share, computed
   from a *different* table. It reproduces to 4 decimals: 0.6818 vs 0.6783. Gated as
   `kmer_share_cross_check_max_diff: 0.01`.
2. **The contrast is nearly identical in magnitude to the GC-vs-dinucleotide contrast
   (+0.270).** Changing the *model class* matters as much as changing the *negative set* —
   and no published benchmark reports both.

---

# Act VI. The verifier's own failure

This is the most important engineering finding in the file, because the verifier is the
project's headline artifact.

## VI.1 The test nobody had run

The reproduction advisor did what should have been done on day one: **corrupted the inputs and
checked that the verifier failed.**

| corruption (on a scratch copy) | checks that FAILED |
|---|---|
| `rehearsal_binding_gc.composition_auroc −= 0.05` | **0 — 71/71 still PASSED** |
| zero every `variant_scores.delta` **and** every conservation value | **0 — 71/71 still PASSED** |
| `variant_specificity.auroc_matched += 0.30` | **0** |
| `cost_of_matching.auroc_dn += 0.05` | 3 |
| `matched_four_models.splicebert −= 0.20` | 4 |
| `robustness.value × 0.5` | 4 |
| `variant_ladder.auroc × 0.5` | 5 |

## VI.2 Why

The gates read **summary tables** (`robustness.csv`, `variant_ladder.csv`) written by the same
analysis pass, not the per-variant sources. Corrupt an input without re-running analysis and
the outputs are unchanged, so verification passes.

**It is a regression detector, and a real one, for the 12 tables it reads.** It is blind to
the 12 it does not. It detects *"did the output change?"*, not *"is the output right?"*

## VI.3 The related provenance hole

`results/tables/variant_scores.csv` is **not** the SpliceBERT arm — it is the **k-mer** arm
(`rehearsal_variants.py:150`). The SpliceBERT matched and mismatched per-variant deltas exist
only in GCS. Consequences:

- Pooled AUROC recomputable from the shipped repo is **0.5765**, not the published 0.829.
- **No headline variant number can be recomputed from the shipped artifact** — including by
  the author.
- `results/` is gitignored (`.gitignore:11`), so a cloner sees zero tables and zero figures.

## VI.4 What it is not

Everything the verifier *does* read reproduced exactly. An independent recompute, with the
advisor's own code and its own RNG, derived the share formula from first principles and
matched all twelve headline numbers:

| quantity | recomputed | published |
|---|---|---|
| R1 GC share | 0.9483 [0.9214, 0.9738] | 0.948 [0.921, 0.974] |
| R1 dinuc share | 0.6783 [0.6207, 0.7347] | 0.678 [0.621, 0.738] |
| R1 difference | 0.2700 [0.2290, 0.3100] | 0.270 [0.231, 0.309] |
| cost of matching | −0.109498 | −0.1095 |
| datasets falling | 94/94 | 94/94 |
| 1-Mb prevalence, ≥20 path | 0.81393 | 0.8139 |
| model − prevalence | −0.05864 | −0.0586 |
| model wins / p | 15/44, p=0.00704 | 15/44, p=0.007 |
| four-model ladder | 0.62793 / 0.68746 / 0.70632 / 0.80908 | 0.6279/0.6875/0.7063/0.8091 |
| SpliceBERT share | 0.4137 | 0.414 |
| size rho | 0.14106, p=0.17506 | 0.141, p=0.175 |

Leave-one-out was verified as genuine: within each block each label class holds exactly one
distinct value, so no variant sees its own label. NaN fraction 8/82 = 0.09756, matching the
golden 0.0976 exactly.

## VI.5 A retraction, in our favour

One advisor flagged that pooled ladder n=18,998 against paired panel Σn=32,353 was
unreconciled, and concluded the 0.829 → 0.755 Simpson correction was "not established."

Asked to dig in, it **retracted**. `variant_ladder` calls `drop_duplicates("vid")` — a variant
near several proteins' peaks is one observation, not several. 32,353 rows → 18,723 unique vids;
the residual 275 is the ladder's slightly wider keep set. Fully explained.

And the correction is **understated**. On a fixed panel (all 82): pooled 0.8294 vs per-dataset
mean 0.6645 = **+0.165** inflation. The same contrast measured on the two arms available
per-variant: conservation +0.032, k-mer +0.021. **The matched arm's pooling inflation is 5-8×
larger than either reference**, exactly because it collapses on the small datasets that pooling
down-weights.

Bug 30 is genuine, and larger than first reported.

---

# Act VII. Bugs 32-50

Continuing `55-the-bug-chronicle.md`, which ended at 31.

## The one that changed the science

**Bug 32 — The negative control's donor was systematically the weaker model.**
*Symptom:* R4's specificity gap was real-looking and passed every gate. *Mechanism:* one donor
per target at a fixed offset on an alphabetical manifest; donors came out weaker than targets
(0.802 vs 0.850, p=0.018), and the gap tracked donor training volume at rho=−0.533.
*Cost:* the paper's one novel positive claim, and ~$4 plus a day to re-establish conditionally.
*Fix:* five quality-spanning donors per target. *Class:* **a control that was never checked
for exchangeability on the covariate that drives the outcome.** The arms were matched on
architecture and nothing else.

## Gates that certified the wrong thing

**Bug 33 — A contamination metric that measured donor size.** `shared_frac = |t∩d|/|t|` is
normalised by the target, so a small donor cannot score high. rho=+0.696 with donor set size;
the flagship correlation reversed and died under adjustment. The size-symmetric `jaccard` was
in the same table, unused. *Class:* an asymmetric ratio used as a symmetric similarity.

**Bug 34 — A golden gate asserting the opposite of its claim.**
`floor_overlap_rho_max: 0.45` capped the floor-vs-contamination correlation. The clean-donor
claim *requires* that correlation to exist; the gate passed as long as it was weak. *Class:*
**a gate that passes when the claim fails is worse than no gate.** This is Bug 29's shape,
one layer up.

**Bug 39 — The verifier passes on zeroed inputs.** See Act VI. *Class:* a regression detector
described as a correctness proof.

## Comparisons on non-identical data

**Bug 37 — The trivial baseline and the model were scored on different variants.**
`cloud_analysis.py`: `auroc_prev` used `s.label[ok]` (blocks with >1 variant) while
`auroc_matched` on the very next line used all of `s.label`. **Mean 22.9% of variants dropped,
max 46.8%.** The "smooth decay" across 100 kb / 1 Mb / 10 Mb was partly the evaluated set
changing size (17,934 / 18,762 / 18,994). *Bias, bounded:* conservation is scoreable both ways
and gives 0.8921 vs 0.8904, so −0.0017. *Corrected headline:* **−0.0605, 15/44, p=0.0044** —
slightly stronger than the −0.0586 it replaced. *Class:* two arms, one mask.

**Bug 38 — Negatives inherit the positive's strand, so ~45% are antisense.** See Act IV.
*Class:* a field deliberately dropped upstream, then needed downstream.

**Bug 40 — The shipped variant table is the wrong arm.** `variant_scores.csv` is the k-mer
arm; pooled AUROC recomputable from the repo is 0.5765 against a published 0.829.
*Class:* a filename that describes a category rather than a member.

## Entrypoints and CI

**Bug 35 — The Modal cost guard had never once executed.** `cloud/modal/guard.py` had
`sys.path.insert` on line 3 and `import sys` on line 47, with `import os` absent entirely.
`python cloud/modal/guard.py` → `NameError: name 'sys' is not defined`. It is the cost control
for the platform where ~95% of the money went, and `docs/56` documented its output as a live
capability. The only broken entrypoint of 11. *Fix:* imports moved above the path insert, `os`
added, docstring made the real module docstring. *Class:* a file nothing imports and no test
runs.

**Bug 36 — CI cannot pass, and has never run.** `ci.yml` installs `docker/requirements-cpu.txt`
(no torch) then runs `pytest tests/unit -q` with no ignores. `test_models.py` imports torch
directly; `test_train_folds.py` pulls it transitively via `src/rbp/train/data.py:10`.
Simulated: 2 collection errors, exit 2. **The fix already existed** in
`docker/cloudbuild.cpu.yaml:55-56` and was never copied over. There is no git remote, so the
workflow has never executed anywhere — and its own opening comment says *"'575 tests pass' is
a claim until something other than the author runs them."*

## Documentation and hygiene

**Bug 41 — Two code comments were false statements.** Both claimed the manifest is
pair-rank-sorted; it is alphabetical (82/82 vs 1/82). One asserted as a design guarantee the
exact property the design lacked. Deleted rather than corrected.

**Bug 42 — The composition share was described as a four-model quantity.** It is measured
against the k-mer model only. The repo was honest (f1's y-axis); the summaries were not.
*Class:* a statistic named after its panel instead of its comparator.

**Bug 43 — `README.md` is still at the pre-run commit.** `git log -1 -- README.md` → `e2ba12f`
*"nothing run yet"*. Publishes R4 at 0.829 (retracted), headlines R3 (cut), and 3 of 4 numbers
in its R4 row match no table. Every correction commit updated `docs/` and skipped the front page.

**Bug 44 — The bug chronicle has four different counts of itself.** Docs say "sixteen"
(55:1, 00:13, 52:272), "Twenty-four" (00:16), "the 24 bugs" (AGENT-CONTEXT:151); the file
contains **25** entries numbered 1-18 and 25-31 — 19-24 do not exist. Also "56/56 golden
checks" (actual 71) and "575 tests" (actual 576). *Class:* **Bug 2 recurring inside the
document written to teach Bug 2.**

**Bug 45 — Live billing account ID in 4 tracked files** while `terraform.tfvars.example`
masks it and the gitignored `terraform.tfvars` says *"Gitignored: it names a billing account."*
The repo states the policy and violates it. A test forbids the *less* sensitive project id and
none forbids this.

**Bug 46 — `results/` is gitignored.** `git ls-files results` returns 0. Defensible for a
pipeline, wrong for a portfolio, and fatal for the claim that anyone can check the 71 numbers.

**Bug 47 — `f4_variant_ladder.png` is 83 minutes older than its own table** and shows the
pre-`abs()`-fix coefficients (~0.71/0.44 against 1.124/0.638). `verify.py` checks tables, not
figures.

## Data-layer findings

**Bug 48 — No ClinVar review-status filter.** `clinvar.py` parses CLNSIG/CLNVC/GENEINFO;
`CLNREVSTAT` has **zero** occurrences repo-wide. 0-star and single-submitter assertions are
included. Exact CLNSIG matching does correctly exclude VUS and Conflicting.

**Bug 49 — Region-class asymmetry in negative sampling.** Positives are classified by
`classify()`'s priority order; negatives are drawn "anywhere in the merged pool." Points drawn
from the cds pool classify as cds only 77.3% of the time. End to end, **6.1-8.7% of negatives
carry a label `classify()` would not give them.**

**Bug 50 — A trivial baseline that is not in the ladder.** Benign and pathogenic ClinVar sets
have different mutation spectra: unique-variant transition fraction **0.705 benign vs 0.517
pathogenic**, with C>A +6.3pp, G>T +5.8pp enriched and T>C −7.9pp, A>G −6.6pp depleted — the
canonical splice-dinucleotide signature. A leave-one-out **substitution-type prior**, no
sequence and no model, scores pooled **0.574**, and on datasets with n≥20 it *beats* the
shipped k-mer delta (0.566 vs 0.551). Must be added to the baseline table.

*Good news in the same family:* pathogenic variants sit closer to peaks (inside 62.0% vs
52.3%, p=5.6e-144) and `−peak_distance` alone gives 0.584 — a nuisance baseline the model
**does** beat, 42/44, p=1.9e-11.

## What this batch has in common

Bugs 1-24 were mostly **infrastructure failing loudly**. Bugs 25-31 were mostly **the science
being subtly wrong while every gate passed**.

**Bugs 32-50 are almost entirely a third shape: artifacts that were correct about the thing
they measured and wrong about the thing they were believed to measure.**

- `shared_frac` correctly measured target-normalised overlap; it was believed to measure contamination.
- `floor_overlap_rho_max` correctly bounded a correlation; it was believed to defend a claim that needed the opposite.
- `verify.py` correctly detects output regression; it was believed to prove correctness.
- `variant_scores.csv` correctly holds the k-mer arm; it was believed to hold SpliceBERT's.
- The single-donor control correctly measured a difference between two heads; it was believed to isolate protein identity.
- The 1-Mb baseline correctly scored the variants it could score; it was believed to be paired with the model.

**No amount of testing finds this class, because every component passes its own test.** What
finds it is asking, for each artifact, *"what would this number look like if the thing I
believe were false?"* — and then computing that. Every one above was caught by a placebo, a
corruption test, a differential, or a cross-check against a second table. None was caught by
an assertion.

---

# Act VIII. Every code change

## New files

| file | lines | what it does |
|---|---|---|
| `scripts/donor_draw.py` | 168 | builds the multi-donor manifest; spans donor quality, ≥2 stronger donors per target, jaccard ≤ 0.02 screen |
| `scripts/multidonor_analysis.py` | 230 | three estimators, target-clustered bootstrap, both panels always |
| `scripts/strand_audit.py` | 174 | GTF gene-strand index; measures the antisense fraction and tests whether it predicts the headline contrast |

## Modified files

```
 cloud/modal/guard.py          |  13 +++--
 cloud/modal/modal_variants.py |  49 +++++++++++++++--
 config/golden.yaml            | 124 +++++++++++++++++++++++++++++++++++++-----
 scripts/cloud_analysis.py     |  73 +++++++++++++++++++++++--
 scripts/variant_splicebert.py |  50 +++++++++++++++--
 scripts/verify.py             | 110 ++++++++++++++++++++++++++++---------
 6 files changed, 359 insertions(+), 60 deletions(-)
```

**`cloud/modal/guard.py`** — imports moved above the `sys.path` insert; `import os` added; the
docstring promoted to a real module docstring. Comment records that the file had never run.

**`cloud/modal/modal_variants.py`** — `donor_task` threaded through `_run_one` and `task`;
`n_donor_tasks()`; `multidonor_sweep` entrypoint using `starmap`; `MAX_CONTAINERS` 10 → 20 with
the cost arithmetic written next to it.

**`scripts/variant_splicebert.py`** — `DONORS` constant; `--donor-task` CLI argument;
`stage_cloud` gained a `donor_task` branch; scoring body factored into `_score(...)` so all
three arms share it; the false pair-rank comment replaced with the measured truth.

**`scripts/cloud_analysis.py`** — three changes:
1. `share()` hoisted from a closure inside `if d is not None:` to function scope, because
   section 1b uses the same estimator on a different table and the two are only comparable if
   it is literally the same function. As a closure, 1b would have been a `NameError` whenever
   `cost_of_matching.csv` was missing.
2. New section **1b**: per-model composition shares plus the k-mer-minus-SpliceBERT contrast,
   with the k-mer row documented as a cross-check.
3. Common-mask columns `n_common`, `auroc_matched_common`, `auroc_conservation_common` so the
   trivial baseline and the model are compared on identical variants.

**`scripts/verify.py`** — `verify_donor_overlap` replaced by `verify_multidonor` (15 new
checks). The docstring records *why*: the old function gated the retracted stratification and
passed 71/71 while the claim underneath was confounded. New gates are deliberately two-sided —
the powered stratum must show the effect **and** the full panel must show its absence.

**`config/golden.yaml`** — 25 new numeric assertions across three blocks:
- `r1_cost_of_matching` gained the per-model shares, the contrast, and the cross-check bound.
- `donor_overlap` → `donor_overlap_RETRACTED`, with the full obituary written in comments
  rather than deleted, because deleting it would hide the most instructive failure in the file.
- `r4_multidonor` (new) — asserts the nulls as well as the effect.
- `strand_audit` (new) — asserts the artifact exists **and** that it does not predict the
  contrast, including `contrast_rho_p_min: 0.05` requiring non-significance.

## New tables

`results/tables/` went from 22 to 28 files:
`donor_tasks.tsv`, `multidonor_pairs.csv`, `multidonor_specificity.csv`, `strand_audit.csv`,
`strand_audit_summary.csv`, plus `variant_tasks.tsv` copied local.

## One self-inflicted error while doing this

Inserting the `strand_audit` block mid-file split `r4_paired_specificity`, orphaning
`size_rho_max` into the new block. `verify.py` caught it immediately —
`verify_r4_paired raised: got KeyError, want no exception ('size_rho_max')` and refused to
certify. Fixed by moving the block to end of file. **The gate worked on its own author**,
which is the first time in this project that has been true of a config change.

---

# Act VIII½. Did the analyst rig the pre-registration?

Asked directly, and it is the right question. Pre-registration only constrains the analyst if
the analyst cannot quietly shape the thing being registered against. Four specific exposures,
each tested rather than argued.

## The structural answer, which is necessary but not sufficient

**The criteria failed.** All three advisors who pre-committed ruled **B3** (conditional), not
**B1** (accept). Two wrote explicitly that they were not softening their own bar. A cheat that
produces a partial negative and demotes the claim out of the abstract is a poor cheat.

That is an argument, not evidence. The evidence follows.

## Exposure 1 — attrition. 475 planned, 409 analysed. Were inconvenient pairs dropped?

No, and it decomposes exactly:

| pairs | cause |
|---|---|
| **65** | 13 targets × 5 donors, all **one-label-class datasets** where an AUROC does not exist |
| **1** | the single failed Modal task, `EFTUD2:HepG2 ← U2AF2:HepG2` |

Those 13 are precisely the known 95 → 82 drop, fixed long before this experiment existed. The
one genuine casualty kept 4 of 5 donors, with gap **+0.0821** against a panel mean of +0.0880
— its loss moves nothing.

Worth recording: the first version of this check flagged the attrition as suspicious, because
it compared median power between lost and kept pairs (86 vs 23) without first separating
whole-target loss from pair-level loss. The heuristic was too crude and misfired on a single
failure. **A check that cries wolf is a check that will be ignored**, so the decomposition is
now by target, not by pair.

## Exposure 2 — did "span, don't match" manufacture the result?

This was the largest exposure. Spanning was *chosen* and matching *rejected*, and that choice
was made by the analyst. So the rejected design was tested directly:

| subset | targets | pairs | gap | wins | p |
|---|---|---|---|---|---|
| all pairs, as reported | 44 | 219 | +0.0880 | 39/44 | 4.6e-08 |
| **tight match** \|Δqual\|≤0.05 **and** \|Δlog size\|≤0.15 | 15 | 16 | **+0.0609** | 11/15 | **0.041** |
| loose match \|Δqual\|≤0.10, \|Δlog size\|≤0.25 | 31 | 46 | +0.1070 | 26/31 | 8.4e-07 |
| donor stronger only (Δqual>0) | 44 | 141 | +0.0920 | 39/44 | 1.3e-07 |
| donor *much* stronger (Δqual>0.05) | 36 | 80 | +0.0971 | 31/36 | 6.3e-07 |

The design that was rejected as infeasible still yields a positive, significant effect on the
16 pairs that happen to satisfy it. Smaller, same sign, same conclusion. **The design choice
did not create the answer** — which also means the feasibility argument for spanning was about
precision, not about getting a different result.

## Exposure 3 — did the "≥2 stronger donors" constraint bias it?

Restricted to the 92 of 95 targets where the constraint was satisfiable: **+0.0913,
p=8.8e-08**. No.

## Exposure 4 — crossed dependence, which an advisor raised as B4 and was not fixed

54 distinct donors serve 219 powered pairs, median reuse 2 and max 44. Target-clustered
bootstrapping does not absorb donor-side dependence.

| clustering unit | intercept | 95% CI |
|---|---|---|
| target (as reported) | +0.0784 | [+0.0512, +0.1034] |
| **donor** | +0.0784 | [+0.0442, +0.1021] |
| conservative envelope, widest bound each side | +0.0784 | **[+0.0454, +0.1043]** |

Still excludes zero. And the structure-free version — **one random donor per target, discarding
the multi-donor design entirely, 2000 redraws**:

```
mean gap +0.0880   2.5-97.5% [+0.0649, +0.1095]   fraction of redraws > 0: 1.000
```

Every redraw is positive. The single-donor design's null was not bad luck at the margin; no
random draw reproduces it.

## What the audit turned up that had been missed

Random single draws give **+0.0880**. The *old* single draw gave **+0.0645** — **smaller**, not
larger. Weaker donors should have *inflated* the old gap, so this needs explaining.

The resolution: the old arm applied **no co-binding screen at all**. Co-binding contamination
raises the floor and shrinks the gap. So the old design carried two opposing biases —
systematically weaker donors (inflating) and unscreened co-binding (deflating) — netting
slightly deflationary.

**Consequence: the `jaccard ≤ 0.02` screen is load-bearing in its own right**, not hygiene. It
belongs in the paper's methods as a substantive design element, and the +0.0645 → +0.0880 move
is partly attributable to it rather than to donor quality.

## The limits that cannot be tested away, stated plainly

**The advisors are instances of the same model as the analyst.** Three agreeing is not three
independent minds agreeing. Correlated blind spots are entirely plausible. This is a
structured self-check with fixed goalposts — genuinely better than unstructured
self-assessment, and genuinely not peer review. Nothing in this document should be read as
having survived external review, because it has not.

**The analyst wrote the brief.** What the advisors saw, and how it was framed, was chosen here.
The unflattering numbers were included and the instruction "do not soften your own threshold"
was explicit, but the frame was still set by the party being audited.

**One B4 concern remains unaddressed:** no attenuation correction on donor quality. It is
measured with error, and regression dilution biases the intercept — the exact quantity being
reported — *toward* the donor-weaker-dominated raw mean, i.e. toward re-creating the confound
that was killed. The direction of that bias is unfavourable and it has not been corrected.

**The strand audit used 40 of 95 datasets** with a regressor of sd 0.038.

## The rule this section exists to state

A pre-registered criterion constrains the analyst only to the extent that the analyst cannot
shape the sample, the estimator, or the exclusions after the fact. Those three are what
Exposures 1-4 test, and all three come back clean. What pre-registration does **not**
constrain is who wrote the question and who chose the reviewers, and in this project both were
the analyst. That residual is real, is not measurable from inside, and is the reason external
review is not optional.

---

# Act IX. The state now

| thing | before | after |
|---|---|---|
| golden checks | 71/71 | **80/80** |
| result tables | 22 | 28 |
| unit tests | 576 pass, 6 skip | unchanged |
| bugs catalogued | 31 | **50** |
| money spent this round | — | **~$4** |
| Modal apps running | 0 | 0 |
| GCP VMs / Batch jobs | 0 | 0 |

## The claim set now

**Headline — R1 model-dependence.** Under dinucleotide-matched negatives, composition
reproduces **0.682** [0.625, 0.736] of a k-mer model's AUROC gain over chance, **0.620** of a
CNN's, and only **0.414** [0.368, 0.458] of SpliceBERT's. Contrast **+0.268 [+0.197, +0.338]**,
nearly identical to the GC-vs-dinucleotide contrast (+0.270).

**Supporting nulls.** phyloP 0.892 and a leave-one-out 1-Mb positional-prevalence rule 0.818
both beat the model 0.755 (paired, common mask: **−0.0605**, 15/44, p=0.0044). Real, and
already published in substance — Grimm 2015, Schreiber 2020.

**Conditional — R4.** Powered 44: +0.088, 39/44, p=4.6e-08; intercept +0.078 [+0.051, +0.104].
All 82: +0.012, n.s.; intercept +0.007 [−0.028, +0.040]. Reported with both panels, labelled
power-conditional, not in the abstract. The novel piece is the detection threshold, not the
effect.

**Stated limitations.** ~45% of negatives antisense (tested; does not create the contrast).
No ClinVar review-status filter. 6-8% region-class mislabelling in negatives. A
substitution-type baseline exists that beats the k-mer delta and is not yet in the ladder.

**Venue,** on which all advisors converged: bioRxiv → *NAR Genomics & Bioinformatics* or
*Briefings in Bioinformatics*. Not *Genome Biology*, not *Nature Methods* — the package is a
benchmarking/protocol paper, and framing it as a research article invites a novelty rejection.

## What remains

| # | item | needs compute | effort |
|---|---|---|---|
| 1 | make `verify.py` **recompute** one headline from per-variant data; ship `scores_sb`/`scores_mm` | no | 2h |
| 2 | rewrite `README.md` against the tables | no | 1h |
| 3 | CI: add the two `--ignore` flags that already exist in `cloudbuild.cpu.yaml` | no | 5m |
| 4 | un-gitignore and commit `results/` | no | 5m |
| 5 | scrub the billing ID from 4 files, add a test | no | 10m |
| 6 | regenerate `f4`; build the missing R5 figure | no | 1h |
| 7 | delete dead stage 10, f3, and R3's orphan golden checks | no | 30m |
| 8 | add the substitution-type baseline to the ladder | no | 30m |
| 9 | fix `negatives.py` strand assignment and re-run | **yes, ~$25** | optional |
| 10 | the manuscript | no | — |

Items 1-8 need no compute, no cloud, and no further review. Item 9 is the only remaining
experiment, and it is optional: the limitation is measured and bounded.

## The one thing to say in an interview about this round

Not "I built a control." The control was wrong. The defensible version:

> I built a wrong-protein control, and it was confounded — my control protein was
> systematically the worse-trained model, and the gap I measured tracked that at rho=−0.53. I
> found it by running a placebo: I split my datasets on a variable that had nothing to do with
> my hypothesis, and it reproduced my flagship result better than the real one. So I re-ran it
> with five quality-spanning donors per target. The confound went to −0.03, and the effect
> survived at +0.088 on well-powered datasets but vanished across the full panel — so I report
> it as power-conditional and keep it out of the abstract. It cost four dollars.

That is a person who can falsify their own result. It is worth more than the result would have
been.

---

# Act X. Rounds 8, 9 and 10: the paper gets chosen, and three of my own results break

Three more council rounds happened after Act IX. They are recorded here because two of them
found errors in work I had already reported as finished, and because the paper's identity
changed.

## Round 8 — which paper, decided by argument

The question was no longer "is this right" but "what is the paper". Four candidates went to a
four-member council: **A** an RBP benchmarking note on R1; **B** a methods paper on three ways
a benchmark delta misleads; **C** a negative-results calibration paper; **D** ship the repo and
stop.

**The editor opened with B**, on the grounds that A carried a revision request — the strand
control — that could destroy its own headline, while B was immune to it. I challenged the
condition they attached, a survey of 40-60 papers, because neither I nor the author can honestly
hand-score papers we have not read; they improved it to a census of papers citing Horlacher 2023.

**Then the prior-art referee killed B.** Every one of its three traps has a published owner, and
two of the remedies are better than mine:

| trap | owner |
|---|---|
| AUROC scale compression | Pencina, D'Agostino, Pencina, Janssens & Greenland 2012, *Am J Epidemiol* 176(6):473-481 |
| coefficients not comparable across fits | Karlson, Holm & Breen 2012, *Sociological Methodology* 42:286-313 — with the significance test mine lacks |
| attenuation against a calibrated null | Schuster, Twisk, ter Riet, Heymans & Rijnhart 2021, *BMC Med Res Methodol* 21:136 |

Whalen, Schreiber, Noble & Pollard 2022, *Nat Rev Genet* 23:169-181 already occupies the
"bundle of genomics ML pitfalls" slot. Put to the editor, they reversed: *"B is dead. Concede
it."* Three of four then said **A**.

**My own contribution to that debate was wrong in an instructive way.** I had measured that two
of the three traps *invert the conclusion* on this data and argued that made B publishable. The
correction is right: that establishes **severity**, not **novelty**, and severity was never in
doubt. I was also overclaiming "2 in 3" as a base rate when n = 3, chosen post hoc from my own
retractions.

**The unlock.** The statistician traced the strand bug to `annotation.py:126`, which states
outright that strand is "deliberately dropped", forcing `negatives.py:328` to inherit the
positive's strand. I read both files and confirmed it. The fix is CPU-only and free; the $25
estimate had always belonged to the SpliceBERT arm, which paper A does not use. A test was
**pre-registered and committed before it ran**: sign retained, CI excluding zero, at least 60%
of the point estimate surviving.

## Round 9 — the test passes, and two of my results break

**The pre-registered strand test passed.** But the design is the interesting part. Restricting
to pairs whose negative is unambiguously sense discards 57% of pairs, and a 256-feature model
loses more from that than a 19-feature baseline does. Restriction alone reports −0.0091 and a
naive reading calls all of it strand. A matched random-drop placebo reproduces −0.0032 of the
same shrinkage with no strand involved.

**Then an auditor showed my placebo was not exchangeable.** Sense-kept pairs are more intronic
(0.4335 against 0.4024 dropped) while GC is balanced (0.5314 against 0.5332), so a uniform
random drop is the wrong counterfactual. With the placebo matched on region:

    restriction alone            -0.0091
     ...cost of dropping pairs   -0.0032
     ...locus mix                -0.0024  [-0.0048, -0.0003]
    STRAND-SPECIFIC (stratified) -0.0036  [-0.0071, +0.0001]
    corrected contrast           +0.0342  90.6% surviving

**A retraction inside a control.** The earlier section said "the artifact is REAL, its interval
excludes zero" on the unstratified −0.0059. Matched on region the interval includes zero. The
artifact is small and *not* distinguishable from zero; only a bound is claimed now. What *is*
distinguishable from zero is the locus-mix component, which is exactly why the plain placebo
was invalid. Two independent implementations agree: −0.0036 mine, −0.0039 theirs.

**And an auditor broke the verifier.** `verify.py` mentioned `rehearsal_binding_gc.csv` and
`rehearsal_binding_dinuc.csv` **zero times**. R1 was gated on `cost_of_matching.csv`, which is a
join of those two, with nothing asserting the copies agree. Permuting the dinucleotide table
passed **166/166** while turning "larger in 88/94" into 67/94, "94/94 fall" into 80/94, and the
Wilcoxon p from 3.8e-17 into 1.5e-12 — the last of which violates this project's *own* threshold,
in a check that never ran on the corrupted table. Eight cross-table assertions at exact equality
now close it; the attack fails at max|diff| 0.282.

**A third finding, cheap and decisive.** The manuscript called the model a **5-mer** in four
places. Both rehearsal tables record **k = 4** on all 189 rows; the "5" came from
`config/params.yaml` `cv: k: 5`, which is the cross-validation fold count. It survived 150 gated
assertions because every one of them checks a *value* and none checked *what produced it*.

## Round 10 — additions, and the cut

The last round asked what could still be added. Everything below was run, not proposed:

- **R1d.** The contrast replicates across cell lines at **r = +0.909** [+0.812, +0.972] over 15
  proteins measured in separate experiments with separately drawn negatives — and it replicates
  *better than either arm it is built from* (+0.518, +0.813). The design guarantees the sign in
  each line independently and guarantees nothing about the magnitudes agreeing, so this is the
  direct answer to "the sign is design-implied". Efficiency rises **×1.31**, so the same
  conclusion needs 58% of the labelled windows.
- **R1e.** The contrast **rebuilt from raw sequence** for all 94 datasets in both arms: +0.0397
  against a committed +0.0397, difference **1.2e-06**. With `recompute.py` this is the second
  thing in the repository that *proves* a headline rather than reproducing it. It is positive at
  every k in 3..6, and k=5 minus k=4 is **+0.0001, p = 0.84** — so the 5-mer error would have
  changed nothing.
- **R1f.** CDS-dominant datasets show **+0.0635** against intron-dominant **+0.0316**,
  p = 1.5e-05, with the mechanism confirmed: intronic sites are more compositional
  (0.6656 vs 0.5765, p = 2.7e-08). The limitation is gated beside it — partialling out total
  gain removes the association, so region indexes *how much* signal exists rather than acting
  independently.
- **The cut.** R3, R4/R4b/R4c, R4d, R5, R6 and the composition-share framing are gone, replaced
  by one table saying what was rejected and why. Manuscript orphans fell **45 → 2**.

**Three defects in my own R4d were fixed even though R4d was being cut**, because the script
ships: the null had used a normal covariate when phyloP is skewed (+1.04), moving it from −68%
to −87%; the script *asserted* that a closed-form approximation "must agree" with the simulation
when that approximation is anti-conservative at this covariate strength and would have fired on
the right answer; and the seed count was 12 against a seed-to-seed sd of 0.054.

## What Act X should be remembered for

Every serious finding in these three rounds was in work I had already reported as complete: a
placebo that was not exchangeable, a verifier that never opened the tables it was certifying, a
model misnamed in four places, and a null built on the wrong distribution. None was found by
thinking harder about the write-up. All four were found by someone running the code against it.

The rule that follows is narrower than "be careful": **a control is not evidence until something
has tried to break the control.** A restriction needs a placebo; a placebo needs stratifying on
whatever the restriction correlates with; a duplicated number needs an assertion that the copies
match; a simulated null needs the real covariate's distribution. Each of those was obvious in
hindsight and none was obvious in advance.
