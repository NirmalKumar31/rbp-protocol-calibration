# 63. The abstract, the title, and what the paper actually argues

**Purpose.** The author is writing the manuscript by hand. This file fixes what it must say, so
the drafting does not have to re-derive any of it. Every number here is gated in
`config/golden.yaml` and reproducible offline with
`python scripts/verify.py --local results/tables` (241/241).

**Register.** Calibration, not bias. The two protocols do not estimate a common quantity, so
nothing may be described as "hidden", "concealed", "inflated" or "proper". Those words
presuppose a true value the paper does not have, and a referee who reads the abstract and then
the limitations will catch the inconsistency.

---

## Working title

> What a benchmark AUROC measures depends on how its negatives were built: a nested
> decomposition across 94 ENCODE eCLIP datasets

---

## Draft abstract, ~330 words

> Sequence models for RNA-binding proteins are benchmarked against negative windows matched to
> the positives on composition, and the choice of what to match on is usually treated as an
> implementation detail. We show it determines most of what the resulting AUROC means.
>
> Across 94 paired ENCODE eCLIP datasets we measure the *nested contribution* of a 4-mer model
> over a 19-feature mono- and dinucleotide composition baseline fitted on the same data. Under
> GC-matched negatives, the protocol in general use, that contribution is **+0.0265 AUROC**.
> Under dinucleotide-matched negatives it is **+0.0662** (difference **+0.0397**, 95% CI
> [+0.0336, +0.0458], larger in 88 of 94 datasets). The same protocol change *lowers* apparent
> AUROC by **0.1095** in 94 of 94. The two quantities a benchmark reports therefore move in
> opposite directions, and a reader ranking protocols by headline AUROC ranks them backwards on
> how much of the model's contribution the benchmark exposes.
>
> Because a nested AUROC gain compresses against a higher baseline, part of the contrast is
> arithmetic. Transporting each arm's increment onto the other's baseline under two link
> functions bounds the protocol effect at **+0.0188 to +0.0313**, positive under every choice.
> A pre-registered strand control, comparing restriction to sense-only pairs against a placebo
> matched on region and gene density, attributes **-0.0055** [-0.0089, -0.0022] to a
> strand-annotation artifact, leaving 85% of the effect. The contrast replicates across cell
> lines at **r = 0.91** over 15 proteins assayed twice, holds at every k from 3 to 6, and is
> twice as large for coding-region binders as for intronic ones (**+0.0635** vs **+0.0316**,
> p = 1.5e-05), consistent with intronic sites being compositionally distinctive.
>
> The two protocols do not estimate a common quantity, so neither figure is the true one. What
> follows is practical: a headline AUROC is uninterpretable without the composition baseline
> measured under the same protocol, and dinucleotide matching measures the model's contribution
> **1.31x** more precisely, reaching the same confidence on roughly 60% of the labelled windows
> for the median dataset. All results are from a single model class; whether they extend to
> deeper architectures is untested.

## Four things about that draft, so they are not lost

1. **The inversion is the hook**, not "matching hides contribution". Apparent AUROC falls in
   94/94 while the nested contribution rises in 88/94. That statement needs no claim that either
   arm is correct, which is exactly what lets it survive the no-common-estimand objection. Note
   in the text that the 94/94 direction is design-guaranteed and all the content is in 88/94, or
   a referee will call the whole thing tautological.
2. **The scope limitation belongs in the abstract**, deliberately. A referee raises it otherwise
   and conceding it costs nothing.
3. **The protocol effect is a RANGE**, +0.0188 to +0.0313. Quoting +0.0313 alone is
   question-begging: the transplant runs both ways and under two links, and the choice moves the
   estimate further than any single interval is wide.
4. **Horlacher et al. 2023 belongs in sentence two of the Introduction**, not the abstract. He
   owns the phenomenon across 11 methods and 313 experiments; this paper owns the nested
   decomposition and the controls.

---

## Section-by-section, what must appear

| section | must contain |
|---|---|
| **Introduction** | Horlacher 2023 as the phenomenon's owner; the gap is what the change is made OF |
| **R1** | +0.0265 / +0.0662, contrast +0.0397 [+0.0336, +0.0458] in 88/94; cost -0.1095 in 94/94; size modification rho +0.307 p=0.0026; the 1-of-15 vs 15-of-15 df concession |
| **R1b** | compression factor 1.51x; the four-member transplant family; d' contrast +0.1290; the log-odds reversal -0.3771 AND its fingerprint (+0.520 vs +0.065) |
| **R1c** | 47.4% demonstrably sense; restriction -0.0091, placebo -0.0040, strand -0.0055 [-0.0089, -0.0022]; corrected +0.0322, 85.4% surviving; the arm asymmetry (+0.047, 37/40, p=5e-09) showing the cue is stronger in the SMALLER-gain arm |
| **R1d** | r = +0.909 [+0.812, +0.972] over 15 proteins; efficiency 1.31x; AND the limitation, partial r = +0.332 [-0.116, +0.690] once total gain is removed |
| **R1e** | rebuilt +0.0397 vs committed +0.0397, difference 1.2e-06; positive at every k, 82/94 at all four; k5 - k4 = +0.0001 p=0.84 |
| **R1f** | CDS +0.0635 vs intron +0.0316, p=1.5e-05; mechanism composition-alone 0.6656 vs 0.5765 p=2.7e-08; AND the limitation, partial rho +0.082 p=0.435 |
| **Discussion** | what a benchmark builder does differently on Monday: report the composition baseline under the same protocol, and prefer dinucleotide matching for precision |
| **Limitations** | design-implied sign; one model class; strand inflates absolute AUROCs in both arms; no expression filter; no eCLIP significance threshold; region-class asymmetry |

## What must NOT appear

R3, R4, R4b, R4c, R4d, R5, R6, and the composition-share framing. All cut, all retained in the
repository, all listed in the manuscript's "Not in this paper, and why" table. Do not resurrect
any of them; each was cut for a reason recorded in `docs/62`.
