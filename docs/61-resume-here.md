# Resume here

Updated 2026-08-28. **241/241 verifier checks, 594 tests, 9 figures, 42 tables, 2 manuscript
orphans.** Branch `r1-scale-check-and-corrected-refit` at `438f3a7`, **22 commits, not merged to
`master`**. Total spend ~$38.03. Working tree clean.

Verify offline with:
`/Users/nirmalkumar/Deep\ Learning\ Project/rna-binding-proteins/.venv/bin/python scripts/verify.py --local results/tables`

## THE ANSWER, IF YOU READ NOTHING ELSE

**The science is finished and it survives scrutiny.** Eleven council rounds; the last one was
unanimous **Minor revision, BLOCKER: NONE** from three independent reviewers who ran code.

**What is left is drafting, and only drafting.** `docs/60-the-paper.md` has all the content
(5,473 words, every number gated) but none of the scaffolding: no Abstract, Introduction,
Results-as-such, Discussion, References, or Data availability. `docs/63` fixes what the abstract
must say so the drafting does not have to re-derive anything. Roughly a day of writing, zero
analysis, zero compute.

Also deferred by the author: README rewrite, billing-ID scrub, git remote and push.

## THE PAPER

**Title (working):** *What a benchmark AUROC measures depends on how its negatives were built: a
nested decomposition across 94 ENCODE eCLIP datasets*
**Venue:** *NAR Genomics & Bioinformatics*, Methods. Realistic outcome: Minor to Major revision.

| section | claim | figure |
|---|---|---|
| **R1** | nested contribution **+0.0265** GC-matched vs **+0.0662** dinuc-matched; contrast **+0.0397** [+0.0336, +0.0458], 88/94; apparent AUROC falls **0.1095** in 94/94 | f0, f1 |
| **R1b** | compression factor 1.51x; protocol effect is a RANGE **+0.0188 to +0.0313** across two transplant directions and two links, positive under all four; d' contrast +0.1290; log-odds reversal -0.3771, diagnosed by the fingerprint (total signal +0.520 vs incremental value +0.065) | f8 |
| **R1c** | strand artifact **-0.0055** [-0.0089, -0.0022], **85.4%** survives; pre-registered before it ran; arm asymmetry shows the cue is stronger in the smaller-gain arm, so the contrast is conservative | f3 |
| **R1d** | replicates across cell lines **r = +0.909**; efficiency **1.31x**. LIMITATION GATED: partial r = **+0.332** [-0.116, +0.690] once total gain is removed | f4 |
| **R1e** | rebuilt from raw sequence, difference **1.2e-06**; positive at every k=3..6; k5 - k4 = **+0.0001**, p=0.84 | f6 |
| **R1f** | CDS **+0.0635** vs intron **+0.0316**, p=1.5e-05; mechanism confirmed. LIMITATION GATED: partial rho +0.082, p=0.435 | f7 |
| R2 | methods table only, citing Horlacher 2023 as replication | f2 |

**CUT and not to be resurrected:** R3, R4/R4b/R4c, R4d, R5, R6, the composition-share framing.
Each has a recorded reason in `docs/62`.

## WHAT CAPS IT, SAID PLAINLY

1. **One model class.** Every number is an L2-penalised logistic on 4-mer counts. R1e varies k,
   not architecture. Stated as a limitation and in the abstract.
2. **The sign is design-implied** (composition spans the 15-df simplex the matcher controls).
   Only the magnitude is informative; R1d is the answer to that.
3. **Horlacher et al. 2023** owns the phenomenon; this paper owns the decomposition and controls.
4. **The biology is weak.** R1f does not survive partialling out total signal, and says so.

This is a solid specialist methods paper. It is **not** a *Genome Biology* or *NAR* paper: a
referee there asks what we now know about RNA-protein recognition, and the honest answer is very
little. The single thing that would change that is showing the contrast holds for a CNN or
SpliceBERT, which needs GPU the author does not have.

## THE DOCUMENT MAP

| file | what it is |
|---|---|
| `docs/59` | council rounds 1-10 as narrative, including Act X |
| `docs/60` | the manuscript's content: every claim, number and limitation |
| `docs/61` | this file: state and what is left |
| `docs/62` | the engineering record: every bug, fix and decision, teaching style |
| `docs/63` | the abstract, the title, and section-by-section requirements |

## STANDING CONSTRAINTS

Ask before Modal or GCS, every time. GCP GPU quota is 0. No em-dashes in the paper or slides.
Keep the Mac awake and cool for long runs (`caffeinate -is`, `OMP_NUM_THREADS=1`, `nice`).
Long runs must be resumable: `strand_placebo.py --resume` keys on an explicit `design` field.

## LESSONS THAT COST REAL WORK. DO NOT RELEARN THEM.

1. **A control is not evidence until something has tried to break the control.** Every serious
   finding in the last three rounds was in work already reported as complete, and every one was
   found by someone running code rather than re-reading.
2. **Restriction is not a random drop.** Any "restrict and recompute" control needs a matched
   placebo, and the placebo needs stratifying on whatever the restriction correlates with, which
   must be MEASURED. Region was not enough; gene density mattered because retention requires
   exactly one overlapping gene strand.
3. **Monte Carlo noise inside a bootstrap inflates the interval.** Five placebo seeds left ~16%
   of the between-dataset variance as seed noise and caused a true finding to be withdrawn. At
   20 seeds it came back. Check the seed count before retracting on an interval.
4. **A duplicated number is worth exactly the assertion that the copies match.** Without it a
   permuted source table passed 166/166.
5. **Never compare a subset rebuild against a full-panel published mean.** This panel-mixing
   error appeared three separate times.
6. **Gates written `if value is not None` disappear rather than fail.** `min_domain_checks`
   fixes the category; it has caught its own author three times.
7. **Check that corroboration is independent before calling it corroboration.** Two derivations
   sharing inputs and algebra are one derivation.
8. **When an estimate needs a modelling choice, report the family, not the member you computed
   first.**
9. **Every one of the 241 assertions checks a VALUE; none checks what produced it.** That is how
   a 4-mer was called a "5-mer" in four places and passed 150/150.
10. **A checker must print its own false-negative rate.** `audit_manuscript.py` is 6.3%
    saturated at four decimals, ~44% at three, and blind to percentages entirely.
11. **Resume logic must key on the design**, not on a column two designs share. The first guard
    silently mixed 38 region-only rows with 2 region-by-density rows.
