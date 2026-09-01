# Resume here

Updated 2026-08-31, after a four-referee council and two rounds of rebuttal. **473/473 verifier
checks, 636 tests, 11 figures, 63 tables, 2 manuscript orphans.** Branch
`r1-scale-check-and-corrected-refit`, **not merged to `master`**. Total spend **~$16 of a $40
budget**, all of it the GC-arm deep sweep; everything since has cost nothing. Working tree clean.

Verify offline:
`/Users/nirmalkumar/Deep\ Learning\ Project/rna-binding-proteins/.venv/bin/python scripts/verify.py --local results/tables`

## THE ANSWER, IF YOU READ NOTHING ELSE

**The science is finished. Stop analysing and start writing.** Four referees, two rounds, no
fatal flaw found. Every remaining idea is a robustness check on a claim that already survived
four attacks; the marginal return went negative around R1o.

**What is left is three things, none of them analysis:** write the manuscript prose, cut the
paper down, scrub the billing IDs and push.

## THE CLAIM, in its final form

> The composition baseline your negative set leaves behind is **most** of what determines a
> sequence model's measurable contribution. Holding model, positives, folds and estimator fixed
> and changing only the negatives, the nested contribution of a 4-mer over a 19-feature
> composition baseline measures **+0.0663 / +0.0265 / +0.0122 AUROC** across three protocols --
> a **5.4-fold range** [4.4, 6.6] -- and falls monotonically as the baseline rises. Given the
> baseline, the protocol label adds **1.0%** of variance; given the protocol, the baseline adds
> **11.0%**. **No rescaling recovers a protocol-free quantity**: over eight monotone transforms
> the range never falls below **2.00x** [1.67, 2.46]. Report the composition-only AUROC under
> the same protocol alongside every headline AUROC.

**Title:** *What a sequence model adds over composition is set by the benchmark's headroom: a
three-protocol calibration across 94 ENCODE eCLIP datasets*
**Venue:** *NAR Genomics & Bioinformatics* or *Bioinformatics Advances*. Editor's forecast:
**~64% / ~72% eventual acceptance.** NOT Genome Biology (~2%) -- that needs a claim about
RNA-protein recognition and this paper has none.

## THE SECTIONS

| section | claim | figure |
|---|---|---|
| **R1** | +0.0265 GC vs +0.0663 dinuc, contrast **+0.0398** [+0.0325, +0.0477] protein-clustered, 88/94; apparent AUROC **−0.1095** in 94/94; helps 80/94 (**72/94** at the measured design effect) | f0, f1 |
| R1b | the transplant family. **SUPERSEDED by R1h**, labelled as such | f8 |
| R1c | strand artifact **−0.0055** [−0.0089, −0.0022], 85.4% survives, pre-registered | f3 |
| R1d | cross-cell-line r = **+0.909**; efficiency 1.31x. Limitation gated: partial r +0.332 n.s. -- but the **log multiplier survives** partialling (+0.680 → +0.580, p=0.038, n=13) | f4 |
| R1e | rebuilt from raw sequence 1.2e-06; positive at every k=3..6 | f6 |
| R1f | CDS +0.0635 vs intron +0.0316. **Dies under its own control** (partial rho +0.082) AND is confounded with achieved match quality. Framing reviewer says cut it | f7 |
| R1g | three model classes: k-mer **+0.0398**, CNN **+0.0530**, SpliceBERT **+0.0864**, 94/94. Multiplier **3.08 / 3.51 / 2.38x** -- the ladder REVERSES on the ratio scale | f9 |
| R1h | the protocol-effect decomposition is **not identified for any model**. Grid: k-mer 4/6, CNN 6/6, SpliceBERT 1/6. The odds link reverses the sign and is named | — |
| R1i | 94 datasets are **79 proteins**; within-protein r **+0.924**; intervals x1.05–1.23; no conclusion changes | — |
| R1j | **40.1%** of negatives untranscribed but BALANCED across arms (p=0.64); excess **−0.0043**, 89.7% survives | — |
| R1k | three protocols, **5.4-fold range**; neg2 lowest at 0.53x | **f10** |
| R1l | protocol and baseline are largely confounded. **CORRECTED:** the 0.0056 AUROC / 3-of-282 figure was the THREE-WAY intersection; pairwise gc-vs-neg2 overlaps 0.212 and 130/188 | — |
| **R1m** | **no rescaling reaches protocol independence**; floor **2.00x** [1.67, 2.46], interval excludes 1. THIS EARNS THE TITLE | — |
| **R1n** | the baseline carries MOST of it: protocol label adds **1.0%** overall, baseline **11.0%**. But decomposed by pair it is 0.05% (gc/dn) vs **3.40%** (gc/neg2), and at matched baseline neg2−gc = **−0.0081 [−0.0130, −0.0036]**, so a protocol-family residual DOES exist. Mechanism: the gradient is a property of composition-matched negatives (−0.545, −0.462) and dies for other-RBPs'-sites (−0.122 n.s.) | — |
| R1o | order-3 baseline removes **64%** of the contrast -- but the fold range survives (2.67x → 2.62x) | — |
| **R1p** | **EXTERNAL VALIDATION.** Horlacher's own negatives, their folds, 45 datasets: range **2.38x** vs our 2.50x, gradient replicates (rho −0.372, p=0.012). **BUT the sign does not reverse there**, so R1n's strong form is ours only | — |

## WHAT WAS WITHDRAWN THIS SESSION, and must not return

1. **"The contrast grows with model capacity"** -- reverses on the ratio scale (R1g) and is
   refuted by R1h's specification grid, where the **CNN** survives most.
2. **"The protocol effect rises with capacity"** -- I "verified" this myself using the very
   assumption R1h refutes. The check was circular.
3. **"The k-mer's protocol effect survives"** -- true only under the most favourable of three
   slope estimators.
4. **"The ordering is not monotone in negative-set hardness"** -- FALSE. Ordered by measured
   difficulty all three orderings are perfectly monotone; neg2 is the EASIEST discrimination,
   not the hardest. I had ranked protocols by "bias-awareness", which is subjective.
5. **"The protocol label carries essentially no information"** -- WITHDRAWN. At matched
   baseline neg2 costs −0.0081 [−0.0130, −0.0036] under two designs. The label carries *little*
   for composition-matched pairs (0.05%) and real information for the other-RBPs'-sites family
   (3.40%). R1p's "failure to replicate" was the same residual, same protocol family.
6. **"protocol effect is +0.0188 to +0.0313"** -- superseded; not identified.

## WHAT IS LEFT

1. **Write the manuscript.** `docs/60` has every claim and number, `docs/63` fixes what the
   abstract must say. Missing: Introduction, Results-as-prose, Discussion, References.
2. **Cut it.** Twelve result sections is at least five too many. `docs/65` carries the framing
   reviewer's specific main-text / supplementary / delete split, and flags ~1,000 words of
   "an earlier draft said X" passages that belong in a response letter, not the paper.
3. **Package for submission.** Billing IDs are now scrubbed and gated by a test, and a LICENSE
   exists. Still to do: git remote and push; a Zenodo deposit; Supplementary Table S1 joining
   the 95 datasets to their ENCFF/ENCSR accessions (recoverable from `config/panel_full*.tsv`);
   a `docs/METHODS.md` lifting the negative matcher's hard-coded parameters (`pool_multiple=8`,
   `pool_min=1500`, greedy cKDTree k=40, L1 over dinucleotide COUNTS) into prose; an AI-use
   disclosure sentence (84 of 101 commits carry a Co-Authored-By trailer); `pdf.fonttype = 42`
   before re-exporting figures, which currently embed Type 3; and f0's ClinVar panel, which
   still advertises a retracted analysis.

**Also missing and it blocks drafting:** `docs/60` has NO write-up for R1j, R1m, R1n, R1o or
R1q -- which includes the result that earns the title and the thesis itself. Their numbers live
only in committed tables and golden.yaml comments.

## STANDING CONSTRAINTS

Ask before Modal or GCS, every time; the one lifting came with a hard cap ($30 credit + $10,
2026-08-29). **GCP billing is CLOSED** -- every object in gs://rbp-repro-2026-derived returns
403, and the pipeline runs off a local store (`rbp.utils.localstore`, `scripts/build_store.py`).
No em-dashes in the paper or slides. Keep the Mac awake and cool (`caffeinate -is`,
`OMP_NUM_THREADS=1`, `nice`). Long runs must be resumable.

**Machine limits, learned the hard way:** 16 GB RAM and ~5 GB free disk. A 945-run local CNN
sweep drove 2.6 GB of swap and a 30x runtime variance on identical work; it was abandoned as
poor value and replaced with a bounded 20-dataset replication that answered the question
(backend effect −0.0030 [−0.0099, +0.0037], contains zero). Bound long local jobs.

## LESSONS THAT COST REAL WORK

1. **A control is not evidence until something has tried to break the control.**
2. **Restriction is not a random drop** -- needs a matched placebo, stratified on whatever the
   restriction correlates with, MEASURED.
3. **Monte Carlo noise inside a bootstrap inflates the interval.** Five placebo seeds hid a
   real strand artifact; twenty found it. **It happened again with expression** (R1j): one seed
   said "no artifact", twenty said −0.0043 with the interval clear of zero. Do not accept OR
   retract an artifact estimate on too few placebo seeds.
4. **A duplicated number is worth exactly the assertion that the copies match.**
5. **Never compare a subset rebuild against a full-panel published mean.**
6. **Gates written `if value is not None` disappear rather than fail.** Same for
   `if col in d else 1.0` -- deleting evidence must not assert the strongest claim.
7. **Check that corroboration is independent before calling it corroboration.**
8. **When an estimate needs a modelling choice, report the family, not the member you computed
   first** -- and NAME the members you exclude (the odds link, R1h).
9. **Almost every assertion checks a VALUE, not what produced it.** That is how a 4-mer was
   called a "5-mer" and passed 150/150, and how the R1g evidence table was forged past 314/314.
   **Fixed by anchoring: both arms' 940 per-window score files are now committed (16 MB), so
   every cell's raw AUROC is recomputable inside verify.py.**
10. **A checker must print its own false-negative rate AND on the range that matters.**
    `audit_manuscript.py` computed saturation over [0.5, 1.0] while R1g's headline lives in
    [0, 0.1], where its 3-dp grid is 100% occupied. It missed 7 of 10 planted lies.
11. **Resume logic must key on the design**, not a column two designs share.
12. **Golden keys built with f-strings are invisible** to `test_golden_keys_are_read`, which
    reads verify.py as TEXT. **Spell them out literally.** This caught its own author again on
    2026-08-31.
13. **A shared venv will import the wrong source tree.** A rootdir `conftest.py` settles it.
14. **Compare rungs on the same rows before comparing them at all.**
15. **A summary table is worth exactly the assertion that it is arithmetic on its evidence** --
    and that check is ONE-DIRECTIONAL. It is equally satisfied by evidence reverse-engineered
    from the summary, which is how the forgery passed. Anchor to per-window evidence.
16. **Traceability is not currency.** A retracted number is still traceable to a committed
    table, so `verify.py` passed 284/284 and the orphan audit reported nothing while the
    LIMITATIONS section quoted three withdrawn strand values.
17. **One-sided floors let inflation through.** `>= 55` passed a rewritten 58 → 94.
18. **Run the test suite BEFORE committing, not after.** Committed a failing build on
    2026-08-31.
19. **An ordering by how principled something sounds is not an ordering.** "Bias-aware" is not
    "harder"; measure the difficulty. This produced a false claim that survived into the
    abstract and was caught by an editor, not by any gate.
