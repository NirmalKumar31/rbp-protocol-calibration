# 60. The paper: every claim, every number, every limitation

**Manuscript source of truth.** Written 2026-08-27 after six rounds of adversarial review.
Every number here is in a committed table under `results/tables/` and gated by
`scripts/verify.py` (105/105). If a number in the manuscript disagrees with this file, this
file is wrong and must be fixed — not the manuscript quietly edited.

**Working title:** *What GC-matched negatives and ClinVar AUROC actually measure: trivial-baseline
calibration of RNA-binding-protein models across 94 ENCODE eCLIP datasets*

**Venue:** bioRxiv immediately; then *NAR Genomics & Bioinformatics* (Methods/benchmarking).

---

## The one-sentence claim

> Under negatives matched only on GC content — the protocol in general use — a 5-mer model's
> measured contribution over a mono+dinucleotide baseline is **+0.0265 AUROC**; matching on
> full dinucleotide composition raises it to **+0.0662** (difference **+0.0397** [+0.0334,
> +0.0461], larger in 88/94 datasets) while lowering the model's apparent AUROC by **0.1095**
> in **94/94** datasets. The standard protocol simultaneously inflates reported performance and
> hides two-thirds of the model's real contribution.

---

## R1 — the protocol effect (primary result)

`cost_of_matching.csv`, n = 94 paired datasets. Figure **f1**.

| quantity | GC-matched | dinuc-matched |
|---|---|---|
| composition alone (19 features) | 0.7827 | 0.6280 |
| composition + 5-mer score | 0.7981 | 0.6886 |
| **nested contribution of the score** | **+0.0265** [+0.0212, +0.0325] | **+0.0662** |
| datasets where the score adds significantly | **80/94** | 82/94 |
| datasets where composition ≥ the k-mer model | **29/94** | — |

- cost of proper matching: **−0.1095** AUROC, **94/94 datasets fall**, paired Wilcoxon p < 1e-15
- contrast in nested contribution: **+0.0397** [+0.0334, +0.0461]
- **effect modification by dataset size is real: rho = +0.307, p = 0.0026** on the reported
  nested contrast. The paired design still rules out size *confounding* — both arms use the
  same datasets — but the effect is larger in larger datasets and that must be printed. (An
  earlier draft printed rho = 0.141, p = 0.175, which is this correlation computed on the
  RETRACTED composition-share difference, not on the quantity the paper claims.)

**The circularity caveat, which must appear in the same paragraph.** The composition baseline's
19 features are mononucleotide + dinucleotide counts + entropy. **GC is not one of them** — it
is the spanned combination C+G. Mononucleotide counts are marginals of dinucleotide counts and
contribute no independent degrees of freedom, so the baseline's frequency space is the 16-cell
dinucleotide simplex: **15 df**. GC is one linear functional of it. So the GC arm controls **1
of 15** and the dinucleotide arm controls **15 of 15**. (An earlier draft said 1-of-18 and
18-of-18 while simultaneously conceding mono is spanned by dinuc; both cannot hold.) The *sign* of the +0.0397 contrast is therefore
implied by the design; only its *magnitude* is informative. Say exactly this. Do not lead with
the 94.8% / 67.8% share framing — see the retraction below.

## R2 — the model ladder (methods table, not a result)

`matched_four_models.csv`, n = 95, identical chromosome-level folds. Figure **f2**.

composition **0.6279** < k-mer **0.6875** < CNN **0.7063** < SpliceBERT **0.8091**.
SpliceBERT beats composition on **95/95**.

**Prior art:** Horlacher et al. 2023, *Briefings in Bioinformatics* published the negative-set
effect across 11 RBP methods. This is scaffolding for R1, cited as replication.

## R3 — the evaluation-protocol asymmetry (new, and it costs nothing)

Recomputed from the 950 committed per-example score files, n = 95 datasets. The design is
1:1 matched pairs, and pooled AUROC discards that pairing.

| model | pooled AUROC | matched-pair concordance | pooling penalty | p |
|---|---|---|---|---|
| SpliceBERT | 0.8091 | 0.8297 | **+0.0206** | 1.8e-16 |
| CNN | 0.7063 | 0.7339 | **+0.0276** | 4.8e-17 |
| **differential (CNN − SpliceBERT)** | | | **+0.0069** | **3.6e-05** |

**A measurement-protocol choice penalises architectures unequally.** It is not an accuracy
artefact: degrading SpliceBERT to CNN accuracy under a binormal decomposition, holding its own
variance structure, predicts a penalty of +0.0065 rather than +0.0276, and rho(penalty
differential, accuracy differential) = -0.14.

**But it is a composition result, not an architecture result, and an earlier draft claimed the
opposite.** The driver is the ratio of pair-shared to within-pair score variance (tau/sigma
0.257 for SpliceBERT vs 0.541 for the CNN, paired p=9e-17; rho(penalty, tau/sigma) = +0.78 and
+0.83). Pairs are matched on composition, so the pair-shared component **is** composition and
locale -- and the penalty correlates with each dataset's composition-only AUROC at rho ~ +0.50.
So R3 is a third view of the R1/R2 composition story rather than an independent finding, and
the sentence "with no composition features involved" has been deleted. Note also that the
retracted pooled ClinVar ladder is the same non-collapsibility phenomenon in the other
direction; presenting one as a result and the other as an artefact requires the explicit
argument that pooled AUROC is unbiased for its own estimand and simply answers a different
question.

## R4 — the model's variant signal is orthogonal to conservation and position (primary ClinVar result)

**RESTRUCTURED 2026-08-27 after peer review.** The previous version of this section led with
"conservation and a positional rule beat the model", which is (a) a replication of Grimm 2015
and Schreiber 2020 and (b) exposed to a circularity the AUROC framing cannot survive: phyloP
is conservation, and ClinVar pathogenic assertions lean on PP3 computational evidence supplied
by conservation-derived tools (CADD, REVEL, GERP), while benign assertions lean on allele
frequency. In this study's own data, benign variants have mean phyloP **0.022** against
pathogenic **5.06**. An AUROC horse race between a model and a partial proxy for the labels is
not interpretable.

**The analysis that is interpretable was already computed and was sitting in a supplementary
table.** It asks a different question — not *who wins*, but *does the model add anything
conservation and position do not* — and that question is immune to the leakage, because
leakage into the labels inflates the competitor's AUROC without creating incremental value for
the model.

Firth logistic regression on 18,830 variants, clustered on 1-Mb genomic blocks,
within-dataset standardised |δ|. `variant_coefficients.csv`, `variant_specificity_refit.csv`.

| arm | controls | coefficient | 95% CI |
|---|---|---|---|
| right protein | none | **1.1238** | [1.032, 1.225] |
| right protein | **+ phyloP** | **1.1214** | [1.036, 1.213] |
| right protein | **+ phyloP + 1-Mb positional rule** | **1.1689** | [1.066, 1.268] |
| wrong protein | + phyloP | 0.6364 | [0.549, 0.712] |
| wrong protein | + phyloP + 1-Mb rule | 0.6918 | [0.623, 0.760] |

- **Attenuation of the model's coefficient after conditioning on conservation: +0.21%.** The
  signal the model reads is essentially independent of phyloP.
- Adding the positional rule as a covariate *raises* the coefficient to 1.1689 rather than
  lowering it.
- The right- and wrong-protein intervals are **non-overlapping** under both control sets
  ([1.036, 1.213] vs [0.549, 0.712]).

**Claim:** a fine-tuned RBP model contributes variant-effect information that neither
conservation nor genomic position supplies, and roughly half of that contribution is specific
to the protein it was trained on. This is stated as incremental value, never as a win.

## R4b — and the two trivial baselines, reported as calibration

Paired per-dataset over the 44 datasets with ≥20 pathogenic variants. Figure **f6**.
**Estimator note:** the rows below are per-dataset means; the pooled versions differ and are
in the supplement. Following peer review, one estimator is used throughout this section — a
paper whose R3 *is* the pooled-versus-paired distinction cannot mix them in its own ladder.

| scorer | AUROC (per-dataset mean) |
|---|---|
| phyloP conservation | **0.892** — but see the circularity above; treat as a ceiling, not a competitor |
| 1-Mb leave-one-out positional prevalence | **0.8139** |
| SpliceBERT, right protein | **0.755** |
| SpliceBERT, wrong protein | 0.691 |

- model vs the 1-Mb rule, common variant mask: **−0.0605**, wins **15/44**, p=0.0044
- decay across block sizes (pooled): 100 kb **0.851**, 1 Mb 0.818, 10 Mb 0.733

**The positional rule is not conservation in disguise, and this matters.** Spearman(prevalence,
phyloP) = **+0.339**, and the rule retains AUROC **0.768–0.832 within every one of ten phyloP
deciles** (unstratified 0.825, n=27,364 unique variants). So it is a second, independent
benchmark leak — and unlike phyloP it is **not** exposed to the PP3 circularity, because no
ACMG criterion is "this variant is near other pathogenic variants." Of the two baselines that
outrank the model, this is the one a referee cannot explain away.

**Prior art, cited in the Introduction rather than buried.** Grimm et al. 2015, *Human
Mutation* (type-2 circularity); Schreiber, Singh, Bilmes & Noble 2020, *Genome Biology*
(average-activity baseline); Livesey & Marsh; Notin et al. 2023 (ProteinGym); Marin et al.
2024 (BEND); Tang & Koo; Sasse et al. 2023; Karollus & Gagneur 2023. **This subsection is a
replication in a new data type.** It is calibration for R4, not a finding.

## R4c — superseded framing, retained for the record

`variant_ladder_paired.csv`, `variant_specificity_attacks.csv`, paired over the 44 datasets
with ≥20 pathogenic variants. Figure **f6**.

| scorer | AUROC |
|---|---|
| phyloP conservation | **0.892** — beats the model 40/44, p = 1.1e-09 |
| 1-Mb leave-one-out positional prevalence | **0.818** (100 kb **0.851**; 10 Mb 0.733) |
| SpliceBERT, right protein | **0.755** |
| SpliceBERT, wrong protein | 0.691 |
| dataset identity alone | 0.668 |
| k-mer \|δ\| | 0.552 |

- paired on a **common variant mask**: model **−0.0605** vs the 1-Mb rule, wins **15/44**,
  p = 0.0044. (The baseline cannot score variants alone in their block — mean 20.2%, max 38.9%
  — so both arms are restricted to the same variants. Bias measured with conservation, which
  is scoreable either way: 0.8921 all vs 0.8904 on-subset, **−0.0017**.)
- the decay across block sizes is smooth, so the rule is reading positional structure rather
  than exploiting one lucky binning.

**Prior art, and it must be cited in the Introduction, not buried.** Grimm et al. 2015,
*Human Mutation* ("type-2 circularity": a predictor knowing only gene identity scores high, and
benchmark rankings invert under gene-disjoint evaluation). Schreiber, Singh, Bilmes & Noble
2020, *Genome Biology* (the average-activity baseline beats deep models that appear to learn
cell-type specificity). Also Livesey & Marsh; Notin et al. 2023 (ProteinGym); Marin et al. 2024
(BEND); Tang & Koo; Sasse et al. 2023; Karollus & Gagneur 2023. **This result is a replication
in a new tissue, not a discovery.** The contribution is the specific instrument (a leave-one-out
positional-prevalence rule on RBP variant data) and the calibration in R5.

## R5 — the wrong-protein control, and the detection threshold (conditional)

`multidonor_specificity.csv`, 95 targets × 5 quality-spanning donors, 474/475 tasks, jaccard
≤ 0.02 screen. Figure **f7**.

| estimator | powered (44) | all usable (82) |
|---|---|---|
| mean gap | **+0.0880**, 39/44, p=4.6e-08 | +0.0124, 56/82, p=0.10 |
| intercept at zero donor advantage | **+0.0784** [+0.0514, +0.1038] | **+0.0068** [−0.0278, +0.0399] |
| ↳ donor-clustered instead of target-clustered | +0.0784 [+0.0442, +0.1021] | — |
| ↳ adding power as a covariate | +0.0236 [−0.0120, +0.0593] | +0.0103 |
| donors *stronger* than the target | +0.0920, 39/44, p=1.3e-07 | +0.0090, p=0.20 |

**Report both panels. The specificity gap does not go in the abstract.**

**The transferable finding is the detection threshold.** The wrong-protein floor is flat in
statistical power (rho **+0.091**, p=0.066) while the model's own arm is steep (rho **+0.631**,
p=7.6e-47). Below ~20 pathogenic variants the target's own head sits at **0.559** — near chance
— so the gap is negative (−0.075, 17/38, p=0.025). **Specificity claims below roughly 20
positives are unmeasurable by construction.** No prior art was found for this calibration.
Figure f7b shows the gap is monotone in power, which is what a real effect does and a threshold
artefact does not.

**Two corrections from peer review, both conceded.** First, the plateau is at **≥45** pathogenic
variants (mean gap 0.108, n=31), not at 20, where the curve is still mid-slope (0.0645). So the
reported +0.0784 is neither a null nor the plateau value; it is an average over whatever power
distribution the 44 happen to have, and the paper's own row — intercept +0.0236 [−0.0120,
+0.0593] once power enters as a covariate — says so. The continuous gap-versus-power curve with
a confidence band, not a stratum mean, should be the primary estimator. Second, f7's x-axis
**filters datasets**, so low-power points are also weaker eCLIP experiments and the 21 points
are nested subsets with dependent p-values; "crosses p<0.05 at 15" is therefore not a test. The
correct experiment holds the 44 powered datasets fixed and **downsamples** pathogenic variants
to n = 5, 10, 20, 40. It needs no GPU and is the third revision we owe. Finally, "the floor is
flat" overstates: the floor rises from 0.634 (<20 positives) to 0.667 (≥20), roughly 6× less
steeply than the matched arm's 0.559 → 0.755, and rho = +0.091, p = 0.066 with CI [−0.01, +0.19]
is a precise near-zero rather than a demonstration of flatness.

**Second methodological finding:** a wrong-protein control must screen its donors on **model
capacity**, not only on protein identity. The v1 design used one donor at a fixed manifest
offset; donors came out systematically weaker (binding AUROC 0.802 vs 0.850, p=0.018) and the
measured gap tracked donor training volume at rho **−0.533**. A placebo split on donor *size*
reproduced the published co-binding stratification **better** (+0.1362 vs +0.1210, both 16/17).
Adebayo et al. 2018 and Hooker et al. 2019 own the general control logic; the RBP-specific
requirement and its failure mode are ours.

## R6 — the substitution-spectrum baseline

`substitution_baseline.csv`. Recomputed here; **weaker than an advisor first reported, and the
conservative reading is the one that ships.**

- 27,492 unique SNVs. Transition fraction: benign **0.705**, pathogenic **0.517** — the
  canonical splice-dinucleotide signature.
- A leave-one-out prior over the 12 substitution types, using no sequence and no model, scores
  pooled AUROC **0.5824**.
- Per dataset (84 datasets with ≥20 pathogenic): substitution prior **0.5627** vs the trained
  k-mer \|δ\| **0.5544**. It beats it on 43/84, mean difference +0.0083, **p = 0.57**.

**Claim it as "neither separable nor equivalent".** p = 0.57 is not evidence of equivalence: the
95% CI on the +0.0083 difference is **[-0.0113, +0.0284]**, a margin wider than half the k-mer
model's entire above-chance signal (0.0544). The honest statement is that both scorers sit near
chance (substitution prior 0.5627 [0.548, 0.578]; k-mer 0.5544 [0.540, 0.568]) and **this panel
cannot separate them**. Do not attribute the class difference to a canonical GT/AG splice
signature: 42.6% of pathogenic variants are C>T/G>A transitions, and the actual class
difference is transversion enrichment (C>A + G>T: 20.3% pathogenic vs 8.3% benign).

---

## RETRACTED, on the record

**The "model-dependent composition share."** Composition reproduces 0.682 of a k-mer model's
gain, 0.620 of a CNN's, 0.414 of SpliceBERT's, contrast +0.268. It is an **algebraic identity**:
`share_m = C / gain_m` where `C = mean(composition_auroc) − 0.5 = 0.127853` is identical across
models, so `share_kmer / share_SB == gain_SB / gain_kmer` exactly — **1.648166 both ways to 6
decimal places**. It is a monotone rescaling of R2, which is prior art, and SpliceBERT beats the
k-mer model on 95/95 datasets, so the contrast excluded zero with probability 1. `verify.py`
now gates the identity rather than a confidence interval. It appears in the manuscript as one
footnote warning against share-as-comparator estimators.

**The pooled ClinVar ladder (matched 0.829).** Simpson-inflated: on a fixed 82-dataset panel
the pooling inflation is **+0.165**, against +0.032 for conservation and +0.021 for the k-mer
model — the deep arm is the one pooling flatters most. Supplement only, labelled as such.

**R3-locality (ISM positional concentration).** Cut: architecture-confounded by construction.

**The "clean donor" stratification (+0.121, 16/17).** The stratifier `shared_frac = |t∩d|/|t|`
is normalised by the target, correlates with donor variant-set size at rho **+0.696**, and the
flagship association reverses and dies under adjustment (+0.299 p=0.049 → −0.240 p=0.12).

---

## LIMITATIONS — manuscript-ready

**1. Negative strand assignment, and the primary result has no strand control.** Sampled
negatives inherit the *positive's* strand while their location is chosen to match composition,
so strand is effectively a coin flip. Audited on 40 of 95 datasets: **55.2%** of negatives carry
the strand of the gene they sit in — but that is the fraction among *unambiguous* windows only.
A further **14.0%** lie in genes transcribed on both strands and are excluded from the
denominator, so the demonstrably-sense fraction is **47.4%**. 0.0% are intergenic. Roughly half
the negative set is therefore antisense sequence no transcript produces, while every positive
is true sense RNA. This inflates all absolute AUROCs in R1.

**An earlier draft defended this with the wrong statistic and the defence is withdrawn.** The
numbers offered (contrast +0.2643 → +0.2787 on sense-only negatives; +0.2485 vs +0.2644 across
panel halves) are the *composition-share contrast* — the quantity retracted below as an
algebraic identity — and one of them appeared in no committed table. R1's reported contrast is
the **+0.0397 nested gain**, and no strand-restricted recomputation of it exists. What survives
from that analysis is narrower but real: the directional cue is detectable and
SpliceBERT-specific (it ranks sense above antisense negatives at 0.576 against the CNN's 0.521,
paired p=3.2e-6), so a pretrained RNA model can read strand. Whether that changes the +0.0397
is **untested**. Regenerating the negatives with strand-correct sampling requires only the
composition and k-mer arms, which are CPU-only, and is the first revision we owe.

*A previous draft also asserted that this convention is shared with the matched-negative
samplers this work benchmarks against. That sentence is deleted: we have no citation for it,
and it is probably false — RNA samplers such as GraphProt draw unbound windows from annotated
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
`src/rbp/data/encode.py` selects by replicate count — the code is right and the documentation
is wrong.

**1d. Conservation is partly a re-read of the labels.** ClinVar pathogenic assertions for
noncoding SNVs rest substantially on PVS1 and PP3, and PP3 is operationalised through
CADD/REVEL/GERP — all conservation-derived; benign assertions lean on BA1/BS1/BP4, i.e. allele
frequency. In this study's own data benign variants have mean phyloP **0.022** (median −0.046,
1.1% above phyloP 7) against pathogenic mean **5.06** (38.4% above phyloP 7). phyloP is
therefore not an independent predictor of these labels but a partial proxy for the evidence
that generated them, and its 0.892 must be read as a **ceiling on this benchmark**, never as a
biological statement. This is the reason R4 is framed as incremental value rather than as a
ranking. The 1-Mb positional rule is **not** subject to this: no ACMG criterion is "near other
pathogenic variants", and the rule holds AUROC 0.768–0.832 within all ten phyloP deciles.

**2. Pretraining exposure.** SpliceBERT is pretrained on human pre-mRNA. Every window in this
study — positive and negative — lies inside an annotated gene (0.0% intergenic), so held-out
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
sharing components gives [0.0215, 0.0998] against [0.0269, 0.0993] iid — the intervals overlap
substantially, so the sharing is real but not materially inflating.

**6. The specificity estimate is not corrected for errors in a covariate.** Donor quality is a
measured binding AUROC with its own sampling error, so regression dilution biases the R5
intercept toward the unadjusted mean. A model-free local estimate on near-neutral donor pairs
(|Δquality| < 0.05, |Δlog size| < 0.25) gives **+0.0774** [+0.0360, +0.1163], 14/18, p=0.0034,
agreeing with the regression intercept to three decimals — so the bias appears small, but it is
not formally bounded.

**7. Review process.** Every claim here has been through six rounds of structured adversarial
critique with pre-registered acceptance criteria, and none of it has been through external peer
review. See `docs/59-the-council-and-the-correction.md`, which records the criteria, the
retractions, and a self-audit of whether those criteria were satisfied by construction.

---

## What the reader should take away

Three things any RBP variant-effect claim must clear before it means anything:

1. **A composition baseline under the negative-set protocol actually used** — under GC matching,
   composition alone reaches 0.783 and the sequence model adds +0.027.
2. **A positional baseline** — a leave-one-out prevalence rule over genomic blocks reaches 0.818
   at 1 Mb and 0.851 at 100 kb, both above a fine-tuned language model at 0.755.
3. **Conservation** — phyloP reaches 0.892 and beats the model on 40 of 44 datasets.

And one design requirement: **a wrong-protein control must match donors on model capacity, and
cannot detect specificity at all below ~20 positive variants.**
