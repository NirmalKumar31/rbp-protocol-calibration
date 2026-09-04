# Execution plan: everything from five reviews, in dependency order

Written 2026-09-04. Sources: four internal referees (statistical, domain, adversarial,
editorial), Gemini, a 64-point external review, `referee_report.md` (M1-M12),
`second_opinion_review.md` (Part 1 audit + N1-N16). Deduplicated against work already done.

**Budget: $20 available. $7 committed to Phase 1a. Hold the rest.**

## THREE CLAIMS VERIFIED AS NON-ISSUES. Do not "fix" these.

1. The nested fit's model-score column **is out-of-fold for every row**, training rows
   included (`_oof_scores` fits on `X[tr]`, whose rows carry OOF scores). Referee M1's
   in-sample worry does not apply. Now stated in Methods.
2. The bias-aware arm reads its donor pool AND target positives from the GC arm's window
   tables. That is why its pair count equals the GC arm's. Disclosed, not a bug.
3. The "reproducibility paradox" (value-assertions do not verify provenance) is already in
   Limitations.

## NUMBERS ALREADY VERIFIED, so no re-derivation needed

* **A1 Table 4.** per-cell mean = 5.417/3.432/2.002/1.476/**1.005**/1.340/1.953 (matches the
  table). Panel-mean = 5.417/3.715/2.548/2.110/1.690/1.488/1.410 (monotone). Per-cell
  **median** = 5.090/3.922/2.608/2.341/1.697/1.365/**1.285**, never reaches 1.
  **62 of 282 cells have (1-c) < 0.15**; min 1-c = 0.0267 (neg2), 0.0335 (gc), 0.1495 (dn).
* **A3.** excess-over-chance span = **18.18**, i.e. it INCREASED the span 3.4-fold.
  Five coordinates reduced, two are identity/affine controls, one increased.
* **A4.** exact binomial 37/40: **9.733e-09** one-sided, **1.947e-08** two-sided. The
  manuscript's 5.0e-9 is some other test and is unnamed.
* **A7.** corr(model-alone AUROC, composition-only AUROC): pooled pearson **+0.894**,
  spearman +0.892. Within arm: dn +0.624/+0.472, gc +0.911/+0.913, neg2 +0.956/+0.949.
* **A5.** contribution range 0.0541 vs standalone gaps 0.019 (CNN) and 0.122 (SpliceBERT).
  Within-dinuc model-class contribution range is 0.1091, i.e. TWICE the protocol range.

---

# PHASE 1 - RECOMPUTATIONS THAT CHANGE NUMBERS

Everything here must finish before any manuscript number is finalised, because each
invalidates tables, golden values and the manuscript audit. Do them together, regenerate
once, re-gate once.

## 1a. Retrain the 20 leaky dinucleotide datasets  [GPU, Modal, ~$7]

The 20: AQR:HepG2 CSTF2T:HepG2 DROSHA:K562 EFTUD2:HepG2 ELAC2:K562 ELAVL1:K562 FKBP4:HepG2
GTF2F1:HepG2 HNRNPM:HepG2 KHSRP:K562 MATR3:K562 NKRF:HepG2 NONO:K562 PCBP2:HepG2 PRPF4:HepG2
PTBP1:HepG2 PTBP1:K562 SLTM:K562 SMNDC1:K562 XPO5:HepG2

**Blockers found, all solvable:**
* There is **no dinucleotide Modal volume**. Only `rbp-gc-store` and `rbp-neg2-store` exist.
  Create `rbp-dinuc-store`, layout mirrors the others: `processed/`, `panel/`, `manifest/`.
* There is **no `sweep_tasks_splicebert_dn.tsv`**. `sweep_tasks_cnn_dn.tsv` exists (476 rows).
  Build the SpliceBERT one the same way.
* `rbp-store/runs/` holds only `gc` and `neg2` - the dn arm predates the local-store era,
  which is *why* these 20 are stale.
* Upload is **91 MB** for the 20 datasets' `dataset.tsv` (whole dn arm is 480 MB).
* Cost: the 20 hold **339,034 of 915,824 rows = 37%** of the arm. The neg2 sweep was 940
  fold-runs / 17.4 GPU-h / $19.10, so 37% ~= 6.4 GPU-h ~= **$7.07**. A run-count scaling
  ($4.06) UNDERSTATES it because these are the panel's largest datasets.
* 200 tasks = 20 datasets x 5 folds x 2 models.

**Procedure, in order:**
1. Back up `data/evidence/scores/{cell}/{protein}/` for the 20, both models, outside the repo.
2. Subset both manifests to the 20 datasets; write `sweep_tasks_splicebert_dn.tsv`.
3. Create the volume; upload the 91 MB + `panel/` + `manifest/`.
4. `RBP_ARM=dinuc modal run cloud/modal/modal_gc_sweep.py::sweep` with an explicit
   dataset list. **The driver currently has `stratified(rows, n)` by COUNT only - add an
   explicit `--only` list, or the resume logic (`done()`, which checks
   `STORE/runs/ARM/.../metrics.json`) will pick the wrong 20.**
5. Smoke-test 2 tasks first. Verify the returned `scores.tsv.gz` `fold` column matches
   `dataset.tsv` before launching the remaining 198.
6. Replace the 20 datasets' committed scores.
7. **Run `scripts/fold_integrity.py` and require `datasets NOT chromosome-grouped, dn arm`
   to be 0, `max chroms/fold` <= 5, and cross-fold neighbours == 0.** Then update
   `golden.yaml`: `grouped_dn: 94`, `leaky_datasets: 0`, `max_chroms_per_score_fold_dn: 5`,
   and add `dn_must_be_fully_grouped: true`.
8. The 74-dataset sensitivity apparatus can then be DELETED from Results and Limitations.

## 1b. Refit the neural arms on the logit scale  [local CPU]  (A24)
The 4-mer enters the nested fit as a log odds, the neural models as probabilities, and
logistic regression is not invariant to a nonlinear transform of a covariate. Measured
effect is small (+0.0008 CNN, +0.0030 SpliceBERT on the GC arm) and it is currently REPORTED
AS A LIMITATION. Fix it instead: refit all neural arms on logit, delete the caveat.

## 1c. Standardise the composition block WITHIN FOLD  [local CPU]  (A25)
**HIGHEST-RISK ITEM ON THIS PLAN.** Currently standardised over the whole dataset. It is
label-free and symmetric so it cannot manufacture a contribution, but it is improper
preprocessing and a certain reviewer comment. Doing it re-derives EVERY composition AUROC,
hence every contribution, hence every table, every golden value and every manuscript number.
Sequence it here or not at all. If deferred, keep the defence paragraph.

## 1d. Regenerate and re-gate, once
Run `run.sh stage 14`-equivalent: every script in `s13b_local_analysis`, then
`audit_manuscript.py`, then `verify.py`. Expect wholesale golden churn from 1c.

---

# PHASE 2 - NEW ANALYSES  [all free, local CPU unless noted]

Ordered so the ones that do not depend on Phase 1 can run in parallel with it.

**Independent of Phase 1 (4-mer only, start immediately):**
* **B5** dinucleotide-**shuffled** fourth arm. Shuffling is what GraphProt/iDeepS/RBPsuite
  actually do and what Tourne indicts. Shuffling preserves composition exactly while
  destroying motif structure, so it should be EASY for a motif model and IMPOSSIBLE for a
  composition model - potentially the one arm where difficulty and contribution move
  together in our own data. ~2 h
* **B6** chromosome-partition sensitivity: 3 alternative partitions meeting the same
  criteria, 4-mer only, report the range of the headline contrast. ~2 h  [candidate for GCP]
* **B8** score the 4-mer standalone on Horlacher's negative sets. **Makes the external test
  like-for-like** - currently their difficulty is read off the composition baseline while
  ours is read off model-alone AUROC, so it never tested the same relation. ~2 h
* **B11** matching-algorithm robustness: greedy vs optimal assignment, pool size 4n/8n/16n,
  on 10-20 datasets. ~2 h
* **B14** show the ~15 nt off-centre signal (asserted in Methods, never shown). ~1 h
* **B15** show the k=3..6 contrast numbers (asserted, never shown). ~30 min
* **B16** summit- vs midpoint-centred window sensitivity. ~2 h
* **B17** transcript-aware region annotation sensitivity. ~2 h
* **B4** bound co-binding label noise. `donor_overlap.csv` already exists. Stratify by
  donor-target peak overlap and show the contribution deficit does not track it. This is the
  most serious untested alternative explanation for the bias-aware arm. ~2 h
* **B10** gene/transcript-clustered CV sensitivity. ~3 h

**Depends on Phase 1 (needs final neural scores):**
* **B1** non-AUROC estimands: Delta deviance / McFadden pseudo-R2, residualised-score AUROC,
  IDI, average-precision increment. **Both long reviews flag that the sweep never leaves the
  ROC**, so "no protocol-free measure" is currently a claim about AUROC only. The penalised
  LRT machinery already exists in `src/rbp/stats.py`. ~3 h
* **B2** promote order-three to a Results subsection: full panel, all three arms, all three
  model classes. Keep the existing caveats (7.16 span is a near-zero-denominator artefact;
  51% of positive mass in 3 of 30 datasets). ~3 h
* **B3** composition-order profile, orders 1-4, per protocol, per model, as the reporting
  object + one figure. Immune to "which order?" and makes the order-three collapse visible.
  ~3 h  [candidate for GCP]
* **B7** 3x3 standalone-AUROC table, all models x all arms. The apparent-AUROC side of the
  thesis currently exists only for the 4-mer. ~1 h
* **B9** rank-normalised-within-fold and fold-averaged AUROC sensitivity (M10). ~1 h
* **B13** interval on r = 0.607 at n = 15. ~15 min

---

# PHASE 3 - MANUSCRIPT FIXES  (30 items, after Phase 1 numbers are final)

## Blocking, numeric/logical
* **A1** Table 4: state aggregation; report the median span; report the (1-c) distribution;
  apply the near-zero-denominator caveat; retract "a monotone rescaling ... does exist".
* **A3** "eight reparameterisations reduced the span" -> five reduced, two are controls, one
  increased it to 18.18.
* **A4** name the test and sidedness behind p = 5.0e-9.
* **A5** Introduction's "nearly three times / about half" - mixed units, and it contradicts
  the Discussion. Make the Intro match.
* **A7** report the +0.894 correlation as the first number in 3.4.
* **A11** pooled Spearman p = 6e-29 treats 282 cells as independent -> clustered interval.
* **A13** multiplicity policy: which analyses confirmatory (adjusted) vs exploratory.
* **A19** Table 1 "within nominal 0.05"; state n per arm (the 1,264-pair gap).
* **A20** "absent" for r = -0.187, p = 0.22 -> "not detectable at this sample size". Text and
  Fig 5b caption.
* **A26** design effect: measured 1.15 / 75-of-94 primary, 1.35 / 72 as sensitivity.
* **A28** explain why the dinuc arm has a LOWER median |dGC| (0.0198) than the GC arm (0.0297).
* **A29** state whether the dinucleotide matcher inherits the >=500 nt peak exclusion.
* **A30** abstract: unfuse the panel-mean/per-dataset sentence; add the order-three clause;
  add "composition beat the 4-mer in 53 of 94 under bias-aware"; state that the
  GC-vs-dinucleotide DIRECTION is design-implied and only the magnitude is informative.

## Blocking, figures
* **A8** Fig 5b groups {gc, dinuc, theirs-n1} vs {neg2, theirs-n2} - the partition 3.4 says
  does not survive. Regroup or caveat.
* **A9** Fig 1b duplicates Fig 3a. Delete one, keep it in Fig 3.
* **A10** Fig 1a/4a use dataset SEs; every headline is protein-clustered. Use those or say why.
* **A18** define Fig 4b's "twelve specifications / six sign tests" grid.

## Blocking, framing and structure
* **A6** title. Candidates: (i) "Measured model contribution is a property of the negative
  set, not the model"; (ii) "How the negatives are built sets what a sequence model appears
  to contribute"; (iii) "Negative-set construction moves measured model contribution
  five-fold and inverts its relation to apparent difficulty".
* **A12** re-head 3.7 to separate proposition A (report the baseline - untested, cannot fail)
  from proposition B (normalise by headroom - pre-specified, failed).
* **A14** Limitations: class ratio 1:1, window size 101 nt, single frozen fold partition.
* **A15** reproducibility reframe: 696 is regression coverage against an expectations file
  generated from the run it validates. Lead with the two end-to-end rebuilds. State that the
  analysis is bit-reproducible while neural score generation is not.
* **A16** verify the Tourne characterisation - "six sampling strategies" and "position weight
  matrix baseline". Two reviewers could confirm neither.
* **A17** extend `audit_manuscript.py` to bare integers; delete the confession sentence.
* **A21** reference style: Horlacher is the most abbreviated entry and is cited ~15 times.
* **A22** resolve the two-family contradiction between 3.4 and 3.6/Discussion.
* **A23** move the partition-defect disclosure from mid-Results to Methods. (Moot if 1a lands.)
* **A27** state why only every second eligible dataset was kept.

---

# PHASE 4 - THE REFRAME  (both long reviews converged on this independently)

* **C1** Lead with: a negative-set protocol is a choice of QUESTION, not a difficulty dial.
* **C2** Consider rebuilding the paper around: **the protocol acts on measured contribution
  entirely through the composition baseline it leaves.** Supported by the caliper-matching
  null (+0.0398 -> -0.0087), by A7's +0.894, and by the within-panel reversal. Crucially it is
  CONSISTENT with Horlacher's benchmark instead of contradicted by it, so it resolves the
  title, the external non-replication and the recommendation's under-specification at once.
* **C3** Package the recommendation as a numbered checklist: negative distribution +
  composition-only AUROC + baseline order.

# PHASE 5 - WRITING AND STRUCTURE
D1 35-40% main-text cut, 3.4 the biggest candidate. D2 sentences over ~40 words.
D3 "placebo" -> "matched negative control". D4 unify neg2/bias-aware/other-RBPs'-sites.
D5 palette consistency + colour-blind check. D6 4 dp is below the per-dataset noise floor.
D7 figure fonts, >=300 dpi vector. D8 `\phantomsection` for back-matter bookmarks.
D9 confirm Supplementary S1 ships.

# PHASE 6 - GCP / MLOPS TRACK  (parallel with everything)
* **E1** `gcloud billing projects link rbp-repro-2026 --billing-account=017994-4FC1D0-8D5176`.
  The account is OPEN and linked to `rbp-composition-2026`; only this project was detached.
  One command makes 25 documented commands work again.
* **E2** re-upload the window store to the live bucket. Enables E3 and un-breaks the docs.
* **E3** run ONE genuinely parallel analysis through `cloud/submit_cpu_sweep.sh` (already
  written, parameterised, documents the 12-vCPU global cap). Use **B6** or **B3** - the
  summary scripts are single-pass and distributing them would be theatre.
* **E4** write the provider split up as a MEASURED decision: T4/L4 quota 0 everywhere, V100
  in 5 regions at 1 preemptible each, `CPUS_ALL_REGIONS` 12 and unraisable, Vertex AI 1 vCPU
  per region so no GPU machine fits -> therefore GPU to Modal at $1.10/GPU-h, CPU-bound
  stages to Batch. This is the strongest MLOps artefact in the repo and it is buried.
* **E5** `terraform plan` clean from an empty project.

# DEFERRED, with reasons
* **F1** anchored negative set (eCLIP SMInput, or RBNS as an orthogonal criterion). 2-3 weeks.
  The only item that would let the paper say WHICH protocol is closer to the truth instead of
  only that they differ. Genome Biology-shaped.
* **F2** a published RBP method on the ladder (RBPNet/GraphProt). 1-2 weeks + GPU.
* **F3** region-matched arm for CNN and SpliceBERT. **~$19**, does not fit alongside 1a.
* **F4** positive-set intersection as PRIMARY. Resisted: real power loss to remove a 0.3%
  discrepancy that is already quantified (Jaccard median 0.9972). Keep as a sensitivity.
* **F5** survey of negative-set reporting in recent RBP papers.
* Zenodo DOI - the user's action. Placeholder is a commented two-line sentence at the end of
  `manuscript/sections/data-availability.tex`. Use the CONCEPT DOI.

---

# BUDGET REVISION, 2026-09-04: $50 available, not $20

The user raised the cap to **$50** explicitly to buy speed. Revised allocation, and an honest
note on what money can and cannot compress.

## What money DOES buy

| item | cost | effect |
|---|---|---|
| 1a retrain, 20 datasets | ~$7 | unchanged; the blocker all five reviews agree on |
| **Phase 2 fan-out onto Modal CPU** | ~$10-15 | **the real win. See below.** |
| **F3 region-matched arm, CNN + SpliceBERT** | ~$19 | now affordable; closes review item N9 |
| headroom | ~$9 | keep |

**THE FAN-OUT IS THE ONLY REAL TIME PURCHASE.** Phase 2 is ~30 hours and most of that is
WAITING on serial single-threaded sklearn loops over 94 datasets on one laptop, on a disk at
98%. Nearly every Phase 2 analysis is embarrassingly parallel per dataset:

* B1 non-AUROC estimands: 94 x 3 arms x 3 models
* B2 order-three, full panel: 94 x 3 x 3
* B3 order profile: 1,128 fits
* B4 co-binding strata, B6 partitions (282), B9 rank-normalised, B10 gene-clustered,
  B11 matching robustness, B16 summit windows, B17 transcript-aware regions

Modal CPU containers are ~$0.05-0.10/core-hour, so all of Phase 2 fanned out is **$10-15**
and turns hours of serial waiting into minutes of wall clock. **Estimated saving: Phase 2
drops from ~30 h to ~15 h.** It also sidesteps the local disk problem.

Practical route: the GPU driver `cloud/modal/modal_gc_sweep.py` already has the manifest +
resume + return-artefacts pattern. Clone that shape into a CPU-only driver
(`@app.function(cpu=..., no gpu=)`), reusing `stratified()`/`done()` and the read-only volume.
Do NOT bend the GPU driver with flags -- `submit_cpu_sweep.sh` documents why two honest files
beat one flagged file.

## What money does NOT buy

The other ~70 hours are manuscript fixes (30 items), the reframe, the 35-40% length cut, and
verification of every change. Those are working time, not compute. **No amount of money
compresses them.** Revised total: **~70-85 h** against the earlier 80-95 h.

## F1 and F2 are still out, and it is TIME not money

* **F1 anchored negative set** (eCLIP SMInput / RBNS): ~25 h plus new data acquisition onto a
  disk with 4.7 GB free. The reviewer's own words: this makes it "a Genome Biology /
  Nature Methods-shaped claim" -- i.e. the next paper, not this revision.
* **F2 published RBP method** (RBPNet / GraphProt): ~20 h of third-party model integration
  with real failure risk (GraphProt is Perl/C++ era).

Compute for both is now affordable; the weeks of engineering are not. Hold as follow-up.

## Revised deliverable

**A1-E5 + F3 + F4 + F5, in ~70-85 h, for ~$36 of $50.**
F1 and F2 deferred to a follow-up.

## DISK: unresolved blocker

4.7 GB free, 98% used. Needs ~4 GB headroom for B16 window re-extraction, the E2 3 GB
staging, and Phase 1 evidence backups. Awaiting user approval to delete:
`rna-binding-proteins/data/raw/GRCh38.primary_assembly.genome.fa.gz` (0.8 GB, redundant with
the uncompressed copy beside it) and `rbp-repo/` (3.0 GB, superseded by rbp-repro).
**Do not delete without asking again.**
