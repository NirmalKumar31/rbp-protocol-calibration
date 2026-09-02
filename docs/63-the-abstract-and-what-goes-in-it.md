# 63. The abstract, the title, and what the paper actually argues

**Purpose.** The author is writing the manuscript by hand. This file fixes what it must say, so
the drafting does not have to re-derive any of it. Every number here is gated in
`config/golden.yaml` and reproducible offline with
`python scripts/verify.py --local results/tables`.

**REWRITTEN 2026-09-02, and the previous version was dangerous.** It still carried the
withdrawn title, it had no R1m, R1n, R1o, R1p, R1q, R1r or R1s, and its item 2 instructed the
author that "the bias-aware protocol is the hardest and reveals the least -- that sentence must
survive editing". **That sentence is false.** Ordered by *measured* difficulty rather than by
how principled the protocol sounds, all three orderings are perfectly monotone and neg2 is the
EASIEST discrimination of the three. A drafting guide that pins a retracted claim and says it
must survive editing is worse than no guide, which is why this file is now regenerated whenever
a claim changes.

**Register.** Calibration, not bias. The three protocols do not estimate a common quantity, so
nothing may be described as "hidden", "concealed", "inflated" or "proper". Those words
presuppose a true value the paper does not have.

---

## Title

> **What a sequence model adds over composition is set by the benchmark's headroom: a
> three-protocol calibration across 94 ENCODE eCLIP datasets**

*Not* "There is no protocol-independent measure of what a sequence model contributes". Two
independent reviewers rejected that for the same reason: a universal negative invites a reader
to go find a counterexample, and three protocols cannot support a universal quantifier. Keep the
negative as the closing line, not the title.

---

## Draft abstract, ~300 words

> Sequence models for RNA-binding proteins are benchmarked against negative windows, and how
> those negatives are chosen is usually treated as an implementation detail. We show it sets
> the scale of the answer.
>
> Across 94 paired ENCODE eCLIP datasets we hold the model, the positives, the
> chromosome-blocked folds and the estimator fixed, and vary only the negative-set protocol. We
> measure the *nested contribution* of a 4-mer model over a 19-feature mono- and dinucleotide
> composition baseline fitted on the same rows. Under GC-matched negatives, the protocol in
> general use, that contribution is **+0.0265 AUROC**; under dinucleotide-matched negatives
> **+0.0663**; under a bias-aware protocol drawing negatives from other RBPs' binding sites in
> the same cell line, **+0.0122**. The same model on the same positives therefore spans a
> **5.4-fold range** [4.4, 6.6], and it falls monotonically as the composition baseline rises.
>
> **The baseline, not the protocol label, is what the range tracks.** Given the baseline,
> knowing which protocol produced it adds **1.0%** of variance; given the protocol, knowing the
> baseline adds **11.0%**. In the 27 of 94 datasets where the bias-aware protocol *lowers* the
> baseline, its deficit reverses (−0.0212 to **+0.0028**), and matched on baseline the headline
> contrast falls from +0.0398 to **−0.0087** [−0.0265, +0.0122]. A protocol-specific residual
> nonetheless survives across protocol *families*: **−0.0081** [−0.0130, −0.0036] at matched
> baseline for the bias-aware arm, whose baseline gradient is itself absent (−0.122, n.s.)
> where composition-matched arms show −0.545 and −0.462.
>
> **No rescaling recovers a protocol-free quantity.** Over eight monotone transforms the range
> never falls below **2.00x** [1.67, 2.46], the minimum being attained by dividing the
> contribution by the baseline's own headroom. On an independent benchmark, using its authors'
> negatives and folds, the range is **2.38x** against our **2.50x**.
>
> Report the composition-only AUROC under the same protocol alongside every headline AUROC.
> Doing so improves cross-protocol rank agreement in 3 of 3 protocol pairs and shrinks
> disagreement in 3 of 3.

## Six things about that draft, so they are not lost

1. **The 5.4-fold range is the phenomenon; the 2.00x floor is the contribution.** Three
   protocols is what makes it a range rather than a contrast, and the third one is what stops a
   referee saying both arms are variants of one flawed design. But R1m is what earns the title:
   without it the paper is "AUROC is compressive near 1", which is a known fact.
2. **Concede the mechanism in the abstract.** The gradient against the composition baseline
   *is* the circularity critique. Stating it first converts it from an objection into the
   finding. A referee who has to discover it will treat it as a flaw.
3. **The ordering IS monotone in measured difficulty.** Do not write "harder negatives reveal
   more" as a paradox, and do not write that the ordering fails to follow stringency. Ordered by
   composition AUROC, apparent AUROC or contribution, all three orderings agree, and neg2 is the
   easiest discrimination. What "harder" buys is headroom, and headroom is set by the protocol.
4. **Do not put the protocol-effect decomposition in the abstract.** R1h shows it is not
   identified for any model. It belongs in the results as a sensitivity band with its own
   failure stated.
5. **Do not write "invariant across model classes".** R1s measures model class at p = 0.023 and
   R1g records SpliceBERT's multiplier as significantly lower. The safe statement is that the
   contrast is present for all three model classes, and that the multiplier is mostly a property
   of the protein.
6. **State the baseline's ORDER wherever a magnitude appears.** R1o shows 87% of the headline
   contrast is one order of composition, so every magnitude in this paper is indexed to a
   mono+dinucleotide baseline. The protocol dependence is not so indexed, and that asymmetry is
   the honest summary of R1m and R1o together.

---

## Section-by-section, what must appear

| section | must contain |
|---|---|
| **Introduction** | Horlacher 2023 *Brief Bioinform* 24(5):bbad307 as the phenomenon's owner, in sentence two; the gap is what the change is made OF |
| **R1** | +0.0265 / +0.0663, contrast +0.0398 [+0.0325, +0.0477] **protein-clustered**, 88/94; apparent AUROC −0.1095 in 94/94; 80/94 helps, **72/94 at the measured design effect**; the 1-of-15 vs 15-of-15 df concession |
| **R1b** | the transplant family, explicitly SUPERSEDED by R1h and labelled as such |
| **R1c** | strand −0.0055 [−0.0089, −0.0022], 85.4% surviving, pre-registered; the arm asymmetry showing the cue is stronger in the SMALLER-gain arm |
| **R1d** | r = +0.909 across cell lines; efficiency 1.31x; AND the limitation, partial r = +0.332, which **repairs on the ratio scale** (+0.580, p = 0.038, n = 13) |
| **R1e** | rebuilt from raw sequence, 1.2e-06; positive at every k = 3..6 |
| **R1f** | CDS +0.0635 vs intron +0.0316; AND that it does not survive partialling total gain, AND that it is confounded with achieved match quality. *A framing reviewer says cut it* |
| **R1g** | k-mer +0.0398, CNN +0.0530, SpliceBERT +0.0864, 94/94; multiplier 3.08x / 3.51x / 2.38x **as a geometric mean over the 77 datasets positive in both arms for all models** -- name the estimator, R1r and R1s use different ones; the ratio-scale REVERSAL and its diagnosis |
| **R1h** | the specification grid: k-mer 4/6, CNN 6/6, SpliceBERT 1/6; no model identified; the odds link named as the member that reverses |
| **R1i** | 94 datasets are 79 proteins; within-protein r = +0.924; intervals x1.05–1.23; no conclusion changes |
| **R1j** | 40.1% untranscribed, BALANCED across arms (+0.000012, p = 0.64); excess −0.0043 [−0.0069, −0.0018]; 89.7% surviving **against this panel's own +0.0414 on 84 datasets, not the published +0.0397 on 94** |
| **R1k** | the three-protocol table; 5.4-fold range; neg2 lowest at 0.53x; the falsified prediction, stated as falsified; the monotone ordering by measured difficulty |
| **R1l** | protocol and baseline confounded by construction; **and that the 0.0056-wide / 3-of-282 common support was the THREE-WAY intersection**, corrected by R1n |
| **R1m** | **the eight transforms and the 2.00x floor** [1.67, 2.46]; Somers' D as the affine control at 5.42x; the excess-normalised coordinate at 18.18x as the normalisation a reader would invent and should not; that the winning coordinate is the paper's own recommendation |
| **R1n** | 1.0% vs 11.0%; the 27/94 natural experiment; matched dn−gc −0.0087; **and the withdrawal** -- pairwise 0.05% (gc/dn) vs 3.40% (gc/neg2), matched neg2−gc −0.0081 [−0.0130, −0.0036], within-arm gradients −0.545 / −0.462 / −0.122 n.s. Two protocol FAMILIES, not three protocols |
| **R1o** | order-3 removes 80/84/88% of the arms and **87% of the contrast**; the fold range does not collapse but **widens, 5.34x → 7.16x**; n = 30; **and the corrected-baseline story**, since the first version measured a different quantity under the paper's name |
| **R1p** | external validation, 2.38x vs 2.50x, gradient rho −0.372 p = 0.012; and that the sign does not reverse there, which R1n's family mechanism now explains rather than concedes |
| **R1q** | rank agreement improves 3/3 and disagreement shrinks 3/3; **only the GC-vs-dinuc rank improvement has an interval clear of zero**, so direction-consistent, not uniformly significant |
| **R1r** | the order-3 baseline absorbs a near-CONSTANT absolute amount from every model (+0.021 gc, +0.054 dn), so the shares are a denominator effect; the result is the per-dataset sign, **k-mer positive in 65/94 against SpliceBERT 94/94**; the protocol contrast survives at order 3 for all three models |
| **R1s** | the multiplier is mostly the protein's: 64.8% against a **29.7% permutation null** for a 79-level factor; cell line noise (p = 0.49); **model class small but real (p = 0.023)** |
| **Discussion** | what a benchmark builder does on Monday: report the composition baseline under the same protocol; never compare contributions across protocols; and if you want to separate model classes, raise the baseline's order rather than the model's capacity |
| **Limitations** | design-implied sign; three architectures not a survey; every magnitude indexed to a mono+dinucleotide baseline (R1o); unseeded initialisation for the deep arms; the CNN rung's CPU-vs-Metal provenance; no expression filter (bounded by R1j); no eCLIP threshold (measured as a non-issue) |

## What must NOT appear

R3, R4, R4b, R4c, R4d, R5, R6, and the composition-share framing. All cut, all retained in the
repository, all listed in the manuscript's "Not in this paper, and why" table.

**And five things that appeared in earlier drafts of this very file or of docs/60. None may
return:**

1. "the protocol effect is +0.0188 to +0.0313" -- superseded, not identified (R1h).
2. "the contrast grows with model capacity" -- reverses on the ratio scale (R1g) and is refuted
   by R1h's specification grid; still lowest for SpliceBERT at order 3 (R1r).
3. **"the ordering is not monotone in negative-set hardness"** -- FALSE. This one was in this
   file, marked as a sentence that must survive editing.
4. "the protocol label carries essentially no information beyond the baseline" -- withdrawn for
   the cross-family comparison (R1n).
5. **"the order-3 collapse shows the 4-mer is uniquely fragile"** -- withdrawn. The absolute
   absorption is the same for every model class (R1r).
