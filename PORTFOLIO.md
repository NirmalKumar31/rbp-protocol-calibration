# How to read this repository

A single-author calibration study, from raw ENCODE files to a submitted preprint: the cloud
pipeline that produced the numbers, the numbers, and the checks that hold them in place.

This file is a reader's guide. `README.md` is the scientific summary; this one says where to
look depending on what you care about, and what each part was hard about.

**No counts are restated here.** Every figure this repository quotes is generated, because a
number typed into a document is wrong by the next commit. `./scripts/preflight.sh` prints the
live ones; `results/tables/release_facts.csv` holds them as data.

## The finding, in one paragraph

Sequence models for RNA-binding proteins are scored against negative windows that somebody has
to construct. Holding the model class, the peak set, the chromosome-blocked folds and the
estimator fixed, and varying only how those negatives are built, a 4-mer model's measured
contribution over a composition baseline moves **4.84-fold** across three protocols. Its
apparent AUROC moves the opposite way. The estimator this literature uses also returns a
positive score for a model whose true contribution is zero by construction; cross-fitting the
score covariate removes at least 95 percent of that. A proposed headroom normalisation was
evaluated against two falsification criteria committed before the external data were scored; it
**failed** on the 45 externally constructed datasets, and the paper reports that failure rather
than dropping it. The study was not formally preregistered, and says so: that commit is a
repository timestamp, not a registry entry.

## Check it yourself, offline, in about a minute

```bash
python -m pip install -e . -c constraints.txt
PYTHONPATH=src python scripts/verify.py --local results/tables
```

That re-derives every published value from committed tables. It is a regression gate, not an
independent reproduction, and `docs/REPRODUCE.md` marks which stages are which.

## Where to look

### If you work on ML systems or infrastructure

Start at **`docs/architecture.md`**. Eight Mermaid diagrams that render in GitHub: the stage
graph and where each stage executes, the two-cloud split as a measured decision, object layout
and the marker discipline, the identity and blast-radius map, the network topology, and one
Batch task end to end.

Then:

| Path | What it shows |
|---|---|
| `run.sh` | The stage graph. Idempotent, resumable, markers written after the payload so an interrupted stage redoes its work |
| `cloud/terraform/` | Custom VPC with no auto subnets, one service account per job, bucket-scoped IAM, budget alert, killswitch |
| `cloud/killswitch/` | Detaches the project's billing account. Reports at startup whether it actually holds the permission, because a dry run proves only the read |
| `cloud/modal/guard.py` | Discovers every sweep app by pattern and sums the bound across all of them |
| `docker/` | Two images. The CPU one ships no torch, which is why the suite has to degrade rather than fail there |
| `scripts/check_image_tree.sh` | Mirrors the Dockerfile's copy set into a temp tree and runs the suite there. The image was unbuildable for weeks before this existed |
| `.github/workflows/ci.yml` | Four jobs, every action pinned by SHA |

The decision worth reading is the GPU split in `docs/architecture.md` section 2: GPU quota was
zero on GCP, so CPU fan-out runs on Batch and GPU work on Modal. The per-model dollar rates were
back-solved from recorded accelerator seconds and then used to price each manifest before
submitting it. `docs/COST.md` separates what was measured from what is forecast.

### If you work on modelling

| Path | What it shows |
|---|---|
| `src/rbp/eval/nested.py` | The estimand: pooled out-of-fold AUROC with the model score added, minus the baseline alone |
| `scripts/cross_fitting.py` | The bias and its removal. Ten extra base fits per dataset to close the outer-fold route |
| `scripts/estimator_floor.py` | The null: a 2-mer whose information the baseline already contains, so the truth is zero |
| `scripts/capacity_ladder.py` | Whether the bias tracks model capacity. It does not; it tracks overlap with the baseline |
| `src/rbp/models/` | 4-mer logistic regression, a DeepBind-style CNN, and a fully fine-tuned SpliceBERT |
| `data/evidence/` | Per-window out-of-fold scores, so every model-class AUROC recomputes from this repository alone |

`scripts/recompute.py` rebuilds published AUROCs from those per-window scores rather than
comparing against a summary, and `scripts/k_sweep.py` rebuilds the headline contrast from raw
sequence. Those two are the only claims here that are re-derived along a different path.

### If you work on analysis and statistics

| Path | What it shows |
|---|---|
| `src/rbp/stats.py` | Protein-clustered bootstrap. Fifteen proteins appear in both cell lines, so resampling datasets would treat correlated rows as independent |
| `src/rbp/eval/delong.py` | DeLong's paired estimator in midrank form, used descriptively because it is not valid against a strictly nested null |
| `scripts/transport_check.py` | All nine train-by-evaluate combinations, decomposed two ways that are different estimands rather than two weightings of one |
| `scripts/estimands.py` | The same comparison on five measures including unbounded deviance, to test whether the result is an artefact of the AUROC ceiling |
| `scripts/scale_check.py` | Eight rescalings, against a simulated equal-means null so the residual span means something |
| `docs/PANELS.md` | Why the dataset counts differ between analyses. The single source; anything disagreeing with it is wrong |

## What is deliberately not here

- **The cross-fitted correction is computed for the k-mer classes only.** The two neural spans
  come from the estimator the paper argues against and are labelled exploratory throughout.
- **An earlier variant-scoring study shares this pipeline.** Its code, its tables and its share
  of the assertions are still here and still pass. `README.md` names every path that belongs to
  it. Nothing in the paper depends on any of it.
- **Two tables no current script reproduces** sit in `results/tables/unattributed/` with a README
  saying so, rather than beside tables that carry a regeneration guarantee.
- **Intermediate window tables are not redistributed**, because they contain genomic sequence.
  They regenerate from the ENCODE accessions in Supplementary Table S1.

## The checks, and what they do not cover

The suite is a regression and provenance gate. Most assertions compare a reported value against
committed evidence from the same run. Two are stronger: published AUROCs rebuilt from per-window
scores, and the headline contrast rebuilt from raw sequence. The manuscript-number trace can find
a value no table supports, but it cannot tell whether a supported value is attached to the right
claim. Scientific validity rests on the design and the limitations, not on the gate count.
