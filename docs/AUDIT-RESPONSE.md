# Response to the external audit of 2026-09-05

An external review of commit `b44e09e` returned 60 numbered items and the verdict "not yet
publication-final". This is the item-by-item response: what was fixed, what was wrong in the
review, and what is deliberately still open.

Its two most useful findings were not on its own P0 list. It measured things nobody had
measured — 193 lint violations, an empty PDF title, a page count of 48 against a documented 28 —
and each of those turned out to be a symptom of one gap: **the checkers had never read the
release documents, and had never read `paper.tex` at all.** The title and abstract were the only
prose in the manuscript no gate had ever traced.

## Resolved

### Scientific

**P0.1, the outer-fold information route.** Correct, and the criticism that a disclosure is not
a measurement was correct. `scripts/cross_fitting.py` reruns the estimator with the route
closed: the covariate on each outer-training row now comes from a base model fitted on the folds
excluding both the outer test fold and the row's own. On the published path the same code
reproduces every published panel mean to `5e-17`, which is what makes the comparison like for
like rather than two implementations.

Two manuscript claims died:

- The Results attributed the estimator's positive floor to *conditioning*. Cross-fitting leaves
  conditioning in place and cuts the floor by 95.7%, 100.7% and 97.6%, landing at `+0.0005`,
  `-0.0001` and `+0.0003` where theory requires exactly zero. The floor is the route.
- The Methods said the bias "can only help the score column". For the 4-mer, closing the route
  *raises* the contribution in all three arms. Withdrawn; the sign depends on the model.

The headline survives and narrows: span 5.42-fold to 4.84-fold.

**P0.2, "truth is zero".** Partly already answered — `results.tex` already gave the
regularisation explanation the review asked for. Superseded in any case: cross-fitting shows the
mechanism is the route, not conditioning.

**P1.1, positive sets differ between arms.** `scripts/common_positives.py` runs the contrast on
the intersection, so the negatives are exactly the only difference.

**P1.6, P1.8, P1.9, P1.12, P1.13.** Five overclaims narrowed to what the design supports:
"invariant to fold size", "unbiased by construction", "could only widen", "independent
benchmark", "specified in advance".

**P1.7.** `nested.py` called DeLong "required" while the Methods cited Demler on its invalidity
for nested models. Both true of different things; the code said only one.

**P1.11.** Every interval is now stated as conditional on six things held fixed, with the
largest of them quantified.

### A finding the review did not make

`scripts/provenance.py` gives every committed table a producing script. Three had none, and one
of those three — `positive_set_overlap.csv` — is quoted in the Discussion. Its `jaccard` column
was `min(n)/max(n)`, exact to `1e-16`: a count ratio labelled a Jaccard index, because it keyed
windows on `id`, which is a per-arm row index. `CSTF2T_pos_0` names a chr12 window in one arm
and a chr9 window on the other strand in the other. Recomputed on genomic coordinate, the median
is unchanged and the minimum moves from 0.9237 to 0.9164.

This is the same failure as the `+0.0397` that motivated `audit_manuscript.py`: a number in the
paper that could not be reproduced and could not fail. That audit closes the
manuscript-to-table link; nothing closed the table-to-script link. Now `run.sh` and CI do.

### Infrastructure

| item | what was wrong |
|---|---|
| P0.5 | Eight stale counts. Root cause: no gate read `README`, `SUBMISSION` or `docs/`, or `paper.tex`. Both audits widened; `release_consistency.py` derives every stated count from the artefact |
| P0.3 | `PROVENANCE.csv` classifies every table raw-reproducible / evidence-recomputable / frozen-only, derived from how `run.sh` invokes it |
| P0.4 | Half wrong: `manifest/study_panel.tsv` is a GCS key, and the study panel *is* committed as `supplementary_table_s1.csv`. Right underneath: the two bias-aware arms had no panel file at all |
| P0.6 | `CITATION.cff`, full package metadata, v1.0.0. The description still named the previous study |
| P1.15 | The hardcoded-ID test watched the *old* project name while the current one reappeared in an argparse default. Now matches the id's shape, exempting comments and docstrings rather than files |
| P1.19–20 | The Modal guard matched one app description while three of four sweeps ran under others, counted one arm of four, and failed **open** on an unreadable CLI |
| P1.21 | `cost.sh` sent gcloud errors to `/dev/null` under `set -uo pipefail`, so an auth failure printed "VMs: none" and "0.00 GB" |
| P1.22 | The killswitch dry run returns before `updateBillingInfo`, so it proved the read while its docstring claimed it proved permissions |
| P2.1 | Ruff 193 to 0 |
| P2.2 | CI ran none of README's claims. It now runs the verifier, the release check, ruff, shell syntax and the LaTeX build |
| P2.4–5 | `constraints.txt`; optional pinned model revisions |
| P2.9 | Five contradictory cost figures were five different quantities. `docs/COST.md` |
| P2.16 | PDF title, author, subject and keywords were empty. Build is now warning-clean |
| P2.19, P2.21 | Per-directory data licences; `results/tables/SCHEMA.md` |

## Wrong in the review

- **P0.4** looked for `manifest/study_panel.tsv` on disk. It is a GCS object key. It also missed
  that the study panel is committed as `supplementary_table_s1.csv`, with the accession and
  inclusion flag per row that the review asked for.
- **P0.2** asked for a caveat that `results.tex:1018` already contained.
- **P2.13** described the abstract as "roughly 255 words" — quoting `SUBMISSION.md`, which was
  stale. It was 616. The review understated this one by trusting a document it was auditing.

## Declined, with reasons

**P2.6, move the narrative comments out of source.** The commentary records why a thing is the
way it is, and in this repository most of it is a bug that cost real time. Moving it to a
changelog puts it where nobody reads it at the moment they need it. It is also several hundred
files of churn against a pipeline whose 990 assertions currently pass.

**P2.18, the `dn` versus `dinuc` naming drift.** Real, and documented in `SCHEMA.md` rather than
fixed: `dn` is a column name that appears in every `check` string in `golden.yaml`, and renaming
it would rewrite the expectations file the checks are made of.

## Still open

Everything here needs new computation or a decision, and none of it is a defect in what is
released.

| | what | why it is open |
|---|---|---|
| P0.1 | Cross-fitting for the CNN and SpliceBERT | four times the GPU sweep, about $76. Stated in the text as not done |
| P1.2 | Multiple negative draws | variability is quantified in the Limitations and not propagated into the intervals |
| P1.3 | Seeded, repeated neural initialisation | needs the GPU sweep rerun |
| P1.4 | Sequence-identity-clustered folds | exact 32-mer sharing and gene-clustered CV are reported; identity clustering is a rebuild |
| P1.5 | Train-fold-only standardisation as primary | measured at `4e-5`; making it primary is a full rerun |
| P1.14 | A 2x2 train-protocol by evaluation-protocol factorial | a separate experiment, and the mechanism claim it would support is not made |
| P1.16–17 | Input checksum manifest and cache lineage hashes | the download code records URL, size and MD5 at fetch time; those records were not committed for the published run, and claiming input integrity now would claim something not held |
| P2.12, P2.15 | 52 pages, 16 table environments | a main/supplement split needs the venue's limits |
| P0.6 | Tag and Zenodo DOI | the author's action, after the content is frozen |

The first four are the next paper rather than a revision of this one: each needs the sweep rerun
that P0.1 needs, and doing them together is one cost rather than four.
