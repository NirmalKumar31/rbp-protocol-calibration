# Resume here

Updated 2026-08-30. **284/284 verifier checks, 615 tests, 10 figures, 44 tables, 2 manuscript
orphans.** Branch `r1-scale-check-and-corrected-refit`, **not merged to `master`**. Total spend
~$54 (~$38 before R1g, ~$16 of Modal credit for the GC-arm deep sweep). Working tree clean.

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
| **R1g** | **the contrast is not a property of the model class and GROWS with capacity**: k-mer **+0.0398**, CNN **+0.0530**, SpliceBERT **+0.0864** [+0.0788, +0.0943] in **94/94**. Paired steps CNN−kmer **+0.0132** [+0.0068, +0.0197], SpliceBERT−CNN **+0.0334**, SpliceBERT−kmer **+0.0466** in 94/94. Protocol effect rises too, +0.0188–0.0314 to +0.0253–0.0543. Refit 4-mer reproduces published R1 to 7.05e-05 | f9 |
| R2 | methods table only, citing Horlacher 2023 as replication | f2 |

**CUT and not to be resurrected:** R3, R4/R4b/R4c, R4d, R5, R6, the composition-share framing.
Each has a recorded reason in `docs/62`.

## WHAT CAPS IT, SAID PLAINLY

1. **Three architectures, not a survey.** RETIRED as the headline limitation by R1g, which
   measures the contrast for a 4-mer logistic, a 7,089-parameter CNN and a 19.7M-parameter
   fine-tuned SpliceBERT and finds it GROWS across all three. What is left is narrower:
   nothing above 20M parameters was tested, so RNA-FM and RNA-MSM (100M, LoRA) are unmeasured,
   and the CNN step is the weak rung (58/94). Do not overstate this in the other direction.
2. **The sign is design-implied** (composition spans the 15-df simplex the matcher controls).
   Only the magnitude is informative; R1d is the answer to that.
3. **Horlacher et al. 2023** owns the phenomenon; this paper owns the decomposition and controls.
4. **The biology is weak.** R1f does not survive partialling out total signal, and says so.

This is a solid specialist methods paper, and R1g moved it up within that band rather than out
of it. The single thing this file previously named as the fix -- showing the contrast holds for
a CNN or SpliceBERT -- **has now been done**, on 2026-08-29, for about $16 of Modal credit plus
2.3 hours of local MPS. It came back stronger than the prediction: the effect is 2.2x larger for
SpliceBERT than for the 4-mer and grows monotonically with capacity.

What that does NOT change: a *Genome Biology* or *NAR* referee still asks what we now know about
RNA-protein recognition, and the honest answer is still very little. R1g makes the
benchmarking claim general across model classes; it does not turn a benchmarking paper into a
biology paper. Aim remains *NAR Genomics & Bioinformatics* or *Bioinformatics Advances*, with a
better shot at the top of that band.

## THE DOCUMENT MAP

| file | what it is |
|---|---|
| `docs/59` | council rounds 1-10 as narrative, including Act X |
| `docs/60` | the manuscript's content: every claim, number and limitation |
| `docs/61` | this file: state and what is left |
| `docs/62` | the engineering record: every bug, fix and decision, teaching style |
| `docs/63` | the abstract, the title, and section-by-section requirements |

## STANDING CONSTRAINTS

Ask before Modal or GCS, every time; the one lifting of that rule came with a hard cap ($30
free credit + $10 out of pocket, 2026-08-29). GCP GPU quota is 0 and the GCP project's BILLING
ACCOUNT IS NOW CLOSED -- every object in gs://rbp-repro-2026-derived returns 403, so the sweep
reads a local store instead (`rbp.utils.localstore`, `scripts/build_store.py`). No em-dashes in the paper or slides.
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
9. **Almost every one of the 284 assertions checks a VALUE, not what produced it.** That is how
   a 4-mer was called a "5-mer" in four places and passed 150/150. R1g's `accelerator` field is
   the counter-example worth copying: it recorded "cpu" for every non-CUDA device until a run on
   Apple MPS exposed it, so the one column a reader would use to check what hardware produced a
   row was wrong by construction.
10. **A checker must print its own false-negative rate -- AND on the range that matters.**
    `audit_manuscript.py` self-reports ~12.5% saturation at four decimals and ~67% at three,
    so the figure quoted here as "6.3% / 44%" was stale by roughly 2x. Worse, it computes
    saturation over **[0.5, 1.0]**, the AUROC range -- but R1g's entire headline lives in
    **[0, 0.1]**, where the 4-dp grid is ~62% occupied and the 3-dp grid is **100%**. Planting
    ten lies in a copy of docs/60: seven were missed, including the central +0.0864 claim. The
    regex also needs a decimal point, so `470/470`, `94/94`, `7,089` and `85.4%` are
    structurally invisible. A self-reported rate measured on the wrong support is worse than
    none, because it is believed.
11. **Resume logic must key on the design**, not on a column two designs share. The first guard
    silently mixed 38 region-only rows with 2 region-by-density rows.
12. **A shared venv will import the wrong source tree.** Several test modules import `rbp` with
    no path setup; the venv's .pth points at `../rna-binding-proteins/src`, and pytest's
    alphabetical order meant the first such module bound the package for every test after it.
    The suite was partly testing a different checkout. A rootdir `conftest.py` settles it.
13. **Compare rungs on the same rows before comparing them at all.** R1g's composition baseline
    check caught a 7th-decimal mismatch (0.59675077 against 0.59675949) that meant the k-mer and
    SpliceBERT were measured on 22,216 and 22,202 rows. Intersect first, fit second.
14. **A summary table is worth exactly the assertion that it is arithmetic on its evidence.**
    Editing `deep_contrast.csv` alone passed every value check until
    `max_summary_arithmetic_diff` existed. Same bug as the permuted `cost_of_matching.csv`.
