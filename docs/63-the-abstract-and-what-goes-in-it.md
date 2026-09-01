# 63. The abstract, the title, and what the paper actually argues

**Purpose.** The author is writing the manuscript by hand. This file fixes what it must say, so
the drafting does not have to re-derive any of it. Every number here is gated in
`config/golden.yaml` and reproducible offline with
`python scripts/verify.py --local results/tables`.

**REWRITTEN 2026-08-31.** The previous version of this file was built around "dinucleotide
matching reveals 2.5x more of a model's contribution than GC matching, and it grows with model
capacity." Both halves are gone. R1h showed the decomposition that supported the first is not
identified; R1k showed the second is false in the only direction that could have tested it. What
replaced them is stronger and is below.

**Register.** Calibration, not bias. The three protocols do not estimate a common quantity, so
nothing may be described as "hidden", "concealed", "inflated" or "proper". Those words
presuppose a true value the paper does not have.

---

## Working title

> There is no protocol-independent measure of what a sequence model contributes: a nested
> decomposition under three negative-set protocols across 94 ENCODE eCLIP datasets

---

## Draft abstract, ~300 words

> Sequence models for RNA-binding proteins are benchmarked against negative windows, and how
> those negatives are chosen is usually treated as an implementation detail. We show it
> determines the answer.
>
> Across 94 paired ENCODE eCLIP datasets we hold the model, the positives, the
> chromosome-blocked folds and the estimator fixed, and vary only the negative-set protocol. We
> measure the *nested contribution* of a 4-mer model over a 19-feature mono- and dinucleotide
> composition baseline fitted on the same rows. Under GC-matched negatives, the protocol in
> general use, that contribution is **+0.0265 AUROC**. Under dinucleotide-matched negatives it
> is **+0.0663**. Under a bias-aware protocol that draws negatives from other RBPs' binding
> sites in the same cell line, it is **+0.0122**. The same model on the same positives therefore
> measures a **5.4-fold range**, and the ordering does not follow negative-set stringency: the
> bias-aware protocol yields the least, at **0.53x** the GC-matched arm.
>
> The mechanism is not subtle and we report it rather than defend against it. Pooled across all
> 282 protocol-dataset combinations, the nested contribution tracks the composition baseline at
> **Spearman −0.60**: most of what a protocol does is decide how much room the baseline leaves.
> Attempts to separate that arithmetic from a residual "protocol effect" by transplanting
> discriminability across baselines are not identified — the transplant's own assumption of
> baseline-invariant increments is violated (slope −0.34, p = 3e-07), and its residuals change
> sign with the direction of transport.
>
> Two controls bound the alternatives. A pre-registered strand control attributes **−0.0055**
> [−0.0089, −0.0022] of the GC-versus-dinucleotide contrast to a strand-annotation artifact, and
> an expression control attributes **−0.0043** [−0.0069, −0.0018] to negatives drawn from
> untranscribed loci, which is 40% of them but is balanced across arms (p = 0.64). Roughly 85%
> and 90% of the contrast survives each.
>
> The practical consequence is a negative result with a positive use: a headline AUROC, and
> even a baseline-relative contribution, is a joint property of a model and a benchmark's
> negative-set construction. Report the composition baseline under the same protocol, and do not
> compare contributions measured under different ones.

## Five things about that draft, so they are not lost

1. **The 5.4-fold range is the claim.** Not "dinucleotide matching is better". Three protocols
   is what makes it a range rather than a contrast, and the third one is what stops a referee
   saying both arms are variants of one flawed design.
2. **The non-monotonicity is the guard against the obvious rewrite.** Someone will want to say
   "harder negatives reveal more". The bias-aware protocol is the hardest and reveals the least.
   That sentence must survive editing.
3. **Concede the mechanism in the abstract.** Spearman −0.60 against the composition baseline
   *is* the circularity critique. Stating it first is what converts it from an objection into
   the finding. A referee who has to discover it will treat it as a flaw.
4. **Do not put the protocol-effect decomposition in the abstract.** R1h shows it is not
   identified for any model. It belongs in the results as a sensitivity band with its own
   failure stated.
5. **The deep models are a limitation section, not a headline.** R1g measured the contrast for a
   CNN and a fine-tuned SpliceBERT and it holds for both (+0.0530, +0.0864), which retires "one
   model class". But the capacity *ordering* is withdrawn on two independent grounds, so the
   only safe statement is "the effect is not a property of the model class".

---

## Section-by-section, what must appear

| section | must contain |
|---|---|
| **Introduction** | Horlacher 2023 *Brief Bioinform* 24(5):bbad307 as the phenomenon's owner, in sentence two; the gap is what the change is made OF |
| **R1** | +0.0265 / +0.0663, contrast +0.0398 [+0.0325, +0.0477] **protein-clustered**, 88/94; apparent AUROC −0.1095 in 94/94; 80/94 helps, **72/94 at the measured design effect**; the 1-of-15 vs 15-of-15 df concession |
| **R1b** | the transplant family, explicitly SUPERSEDED by R1h and labelled as such |
| **R1c** | strand −0.0055 [−0.0089, −0.0022], 85.4% surviving, pre-registered; the arm asymmetry showing the cue is stronger in the SMALLER-gain arm |
| **R1d** | r = +0.909 across cell lines; efficiency 1.31x; AND the limitation, partial r = +0.332, which **repairs on the ratio scale** (+0.580, p = 0.038, n = 13, interval touches zero) |
| **R1e** | rebuilt from raw sequence, 1.2e-06; positive at every k = 3..6 |
| **R1f** | CDS +0.0635 vs intron +0.0316; AND that it does not survive partialling total gain, AND that it is confounded with achieved match quality |
| **R1g** | k-mer +0.0398, CNN +0.0530, SpliceBERT +0.0864, 94/94; multiplier 3.08x / 3.51x / 2.38x; the ratio-scale REVERSAL and its diagnosis |
| **R1h** | the specification grid: k-mer 4/6, CNN 6/6, SpliceBERT 1/6; no model identified; the odds link named as the member that reverses |
| **R1i** | 94 datasets are 79 proteins; within-protein r = +0.924; intervals x1.05–1.23; no conclusion changes |
| **R1j** | 40.1% untranscribed, BALANCED across arms (+0.000012, p = 0.64); excess −0.0043 [−0.0069, −0.0018]; 89.7% surviving |
| **R1k** | the three-protocol table; 5.4-fold range; neg2 lowest at 0.53x; the falsified prediction, stated as falsified; Spearman −0.60 |
| **Discussion** | what a benchmark builder does on Monday: report the composition baseline under the same protocol; never compare contributions across protocols |
| **Limitations** | design-implied sign; three architectures not a survey; unseeded initialisation for the deep arms; the CNN rung's CPU-vs-Metal provenance; no expression filter (bounded by R1j); no eCLIP threshold (measured as a non-issue) |

## What must NOT appear

R3, R4, R4b, R4c, R4d, R5, R6, and the composition-share framing. All cut, all retained in the
repository, all listed in the manuscript's "Not in this paper, and why" table.

**And two things that were in earlier drafts of this very file:** "the protocol effect is
+0.0188 to +0.0313" (superseded by R1h) and "the contrast grows with model capacity" (falsified
on the ratio scale by R1g and refuted by R1h's specification grid). Neither may return.
