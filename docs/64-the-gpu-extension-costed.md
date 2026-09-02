# 64. Lifting the paper past one model class: what it costs, measured

Written 2026-08-29, before anything was spent. Budget available: a fresh Modal account with
$30 credit plus $10 out of pocket, so **$40**.

The limitation in `docs/61` is "one model class. Every number is an L2-penalised logistic on
4-mer counts." This file prices removing it. Every figure below is calibrated against work
this project already paid for, not against a spec sheet.

---

## 1. The position is better than `docs/61` says

Three things were already on disk and had not been noticed.

**The dinuc arm of the deep-model experiment is already done and committed.**
`data/evidence/scores/{cell}/{protein}/{cnn,splicebert}/fold{0..4}/scores.tsv.gz` holds
per-window out-of-fold scores for 95 datasets, both models. Verified: pooling AATF:K562's
five folds reproduces `matched_four_models.csv` at 0.7553366918800132, exact to the last
digit. 475 runs, ~$31 of GPU time, already bought.

**The GC-arm windows are on local disk.** `rna-binding-proteins/data/processed/{cell}/{protein}/dataset.tsv`
is the GC arm for all 94 paired datasets, with `fold` already assigned and `seq_rna` ready.
AATF:K562 is 2086 rows, matching `panel_final_K562_gc.tsv` (1043 pairs) and
`rehearsal_binding_gc.csv` (n=2086). Nothing needs preparing or re-downloading. 244 MB.

**So the only missing piece is GC-arm training.** `cloud/modal/modal_sweep.py` hardcodes
`"ARM": "dinuc"` in `_env()`. That one line is the whole gap.

### The one thing that got worse

**GCP billing is disabled.** `gcloud storage ls gs://rbp-repro-2026-derived/` returns
`HTTPError 403: The billing account for the owning project is disabled in state absent`.
The buckets still enumerate, so the objects are probably intact, but they cannot be read.
The Modal sweep reads datasets from GCS and writes scores back to it, so **that path is dead
until either billing is restored or the data path is replaced**. See item H.

---

## 2. Half the experiment, computed for $0

Ran `gain_over_composition` over the committed dinuc scores, all 94 paired datasets.

| model | nested contribution, dinuc arm | positive in |
|---|---|---|
| composition baseline | AUROC 0.6274 | |
| k-mer (the paper's model) | **+0.0662** | 93/94 |
| CNN | **+0.0860** | 89/94 |
| SpliceBERT | **+0.1754** | **94/94** |

The k-mer row reproduces the published +0.0662 exactly, which validates the pipeline.
SpliceBERT's nested contribution is **2.65x the k-mer's** and positive on every dataset.

**Why this matters for the decision.** It is the dinuc half of the contrast. Only the GC
half needs buying, and the quantity being measured is large and clean rather than marginal.

### What the GC arm would have to beat

Running R1b's own d' transplant onto these numbers, compression alone predicts:

| model | dinuc gain | contrast predicted by compression alone |
|---|---|---|
| k-mer | +0.0662 | +0.0182 (**observed +0.0397**, so protocol effect +0.0215) |
| CNN | +0.0860 | +0.0284 |
| SpliceBERT | +0.1754 | +0.0552 |

The k-mer row is the check: the transplant reproduces `contrast_protocol_reverse = +0.0215`
from `scale_check.csv`. So the machinery is right, and the SpliceBERT row is a real
prediction. Anything the GC arm measures above +0.0552 is protocol effect.

---

## 3. Everything that needs a GPU, priced

Calibration anchor: the dinuc SpliceBERT sweep was 475 runs over 916,672 windows, measured
at **2.33 h wall on 10 concurrent A10G** and recorded at **~$31** (`docs/57`, stage 9).
The GC arm is 913,468 windows, **99.7% of the same work**.

Measured throughput, this project's own batch size and fp32 (no autocast in `trainer.py`):

| | SpliceBERT | CNN |
|---|---|---|
| A10G, from the sweep log | 268 train window-visits/s | |
| **Apple M4 MPS, benchmarked today** | **106 window-visits/s** | **3428 window-visits/s** |

The Mac is only **2.5x slower than an A10G** on SpliceBERT. That makes local a genuine
fallback rather than a joke.

| # | item | needs GPU | Modal cost | local (M4) | verdict |
|---|---|---|---|---|---|
| **A** | **GC-arm SpliceBERT, 94 datasets x 5 folds = 470 runs** | yes | **~$31**, 2.5 h wall | ~55 h compute, realistically 70-90 h | **the whole point** |
| A' | same, size-stratified 40 datasets | yes | ~$14 | ~32 h | cheaper fallback |
| A'' | same, size-stratified 30 datasets | yes | ~$11 | ~20 h | floor |
| **B** | **GC-arm CNN, 470 runs** | optional | ~$3 | **~2-4 h, free** | run it locally |
| **C** | nested contribution, both arms, both models | **no** | $0 | minutes | already proven today |
| **D** | model-ranking inversion across arms | **no** | $0 | free | falls out of A+B |
| **E** | R1d replication across cell lines, for SpliceBERT | **no** | $0 | free | falls out of A |
| F | frozen SpliceBERT probe, both arms | yes | ~$1 | ~2 h | optional, neat |
| G | RNABERT (0.48M), both arms | yes | ~$2 | ~3 h | optional, marginal |
| H | **replace the dead GCS data path with a Modal Volume** | no | $0 | 1-2 h work | **required** |
| X | RNA-FM / RNA-MSM (100M, LoRA) | yes | **~$150 each** | days | **out of budget** |
| X | R3 locality / ISM on the GC arm | yes | ~$0.30 | | R3 is cut, skip |

LoRA does not help X: it backpropagates through the frozen stack, so the FLOPs scale with
parameters regardless. 5x the parameters is 5x the bill.

### The panel-size lever

Cost scales with total windows, and the panel is very skewed (median 6,104, max 32,384).

| subset | share of GC-arm work | Modal cost | expected t on the protocol effect |
|---|---|---|---|
| size-stratified 20 | 25.4% | ~$8 | 5.4 |
| size-stratified 30 | 35.4% | ~$11 | 6.6 |
| size-stratified 40 | 46.4% | ~$14 | 7.7 |
| size-stratified 60 | 67.3% | ~$21 | 9.4 |
| all 94 | 100% | ~$31 | 11.7 |

The t column comes from the k-mer's own per-dataset protocol effect (mean +0.0215, SD
0.0177), which is well determined because the transplant is paired per dataset. Even n=20
clears t=5. **Statistically a subset is fine; the reason to prefer all 94 is that the
headline is n=94 and a referee should not have to ask why the deep models used a different
panel.**

Do **not** subset by taking the smallest datasets, even though it is 5x cheaper. R1's own
size modification (rho +0.307, p=0.0026) says the contrast is larger in bigger datasets, so
smallest-first is the one subset guaranteed to be biased.

Do **not** subsample windows within datasets. The dinuc arm was trained on full datasets;
capping the GC arm would confound arm with training-set size, and re-training the dinuc arm
to match costs more than it saves.

Do **not** change batch size, epochs or precision to go faster. The two arms must be trained
identically or the contrast confounds protocol with hyperparameters.

---

## 4. Recommended plan, $33 of $40

| stage | what | cost | why |
|---|---|---|---|
| 0 | dinuc-arm nested contribution for CNN + SpliceBERT | **$0, done** | already computed above |
| 1 | **H**: Modal Volume data path, no GCS | **$0**, 1-2 h | GCP billing is dead; also makes the repo runnable without GCP, which is better for the MLOps story |
| 2 | **pilot**: GC SpliceBERT on 12 stratified datasets | **~$2** | proves the path end to end, and measures the GC arm's epochs-per-run so stage 3 can be re-priced before committing |
| 3 | **A**: GC SpliceBERT, all 94 | **~$31** | the result |
| 4 | **B**: GC CNN, locally overnight | **$0** | third model class |
| 5 | **C/D/E**: analysis, figures, golden block, verifier assertions | **$0** | |

**Total ~$33**, leaving ~$7 for retries. If the pilot re-prices stage 3 above ~$35, fall back
to the size-stratified 60 (~$21) and say so in the methods.

Do not skip stage 2. Every extrapolation in this project has been wrong: the CNN cloud
penalty was 1.65x, RNABERT's was 4.9x, the GPU speedup was guessed at 100-200x and measured
at 29.6x, and a cost estimate built from three published prices came out 44% high.

---

## 5. What could go wrong

1. **The result could be a null**, i.e. SpliceBERT's contrast lands at or below the +0.0552
   that compression alone predicts. That is not a wasted $31: it is the finding that the
   protocol effect is specific to composition-adjacent models and vanishes for a pretrained
   transformer. It changes the paper's claim from "general" to "specific" and both are
   publishable. **This experiment cannot fail to produce a reportable answer**, which is the
   main reason it is worth buying.
2. **The GC arm's epoch count is unknown.** The GC task is easier (composition baseline
   0.7827 vs 0.6274), so early stopping could fire sooner or later. Budget +/- 20%. Stage 2
   measures it.
3. **Concurrency is the only real cost guard.** `MAX_CONTAINERS = 10` at A10G is $11/h. There
   is no Modal equivalent of the billing killswitch. Leave it at 10.
4. **The GCS objects may not be recoverable at all.** Does not block anything: every input
   needed is on local disk. It only means the 475 dinuc checkpoints are gone, and we do not
   need them, because the scores are committed.
5. **New numbers mean new gates.** Anything added here needs a `golden.yaml` block and
   assertions, and `min_domain_checks` has to rise. Untested gates are assumed broken until
   an attack against them fails (see `docs/62`).

---

## 6. What actually happened, 2026-08-29 to 30

Every estimate in this file was made before anything ran. Here is the outturn, so the next
estimate can be calibrated against a real one rather than another guess.

| | predicted | actual |
|---|---|---|
| GC-arm SpliceBERT, 94 datasets | ~$31, 2.5 h | **~$16.02**, 470/470 runs, 0 failures |
| training GPU-h for that arm | 11.8 (from the pilot) | **11.99** |
| GC-arm CNN | $0 local, 2-4 h | **$0 local, 2.27 h**, 470/470, 0 failures |
| total spend | ~$33 | **~$16.10** of a $40 cap |

**The pilot paid for itself twice over.** It re-priced the full run from $31 to $16 -- the GC
arm turned out to run about 2x faster per pair than the dinucleotide arm did, on the same
epochs (7.56 against 7.68), so the difference is hardware and not science -- and its 59 runs
were kept and skipped by the full sweep, so the checkpoint cost nothing. The linear fit of
seconds against pairs across the 12 pilot datasets had r = 0.9867, which is why extrapolating
from 10% of the work was safe.

**The estimate that was wrong was the pessimistic one**, which is the direction to be wrong in.
The $31 anchor came from the dinucleotide arm's own bill and was honest; the hardware moved
underneath it.

### The result, against the risk stated in section 5

Section 5 said the experiment could come back null and that this would still be publishable.
It came back the other way: the contrast is **+0.0398 (k-mer), +0.0530 (CNN), +0.0864
(SpliceBERT)**, growing monotonically with capacity, positive in 94/94 for SpliceBERT, and the
protocol effect rises with it. The compression-only prediction in section 2 was +0.0552 for
SpliceBERT against an observed contrast of +0.0864, so the transplant's forecast was in the
right place and the excess is protocol.

### What was not foreseen

1. **The GCP billing account was closed, not merely unfunded.** Item H was written as "replace
   the data path" and turned out to be the whole prerequisite. `rbp.utils.localstore` plus
   `scripts/build_store.py` now make the sweep runnable with no cloud storage at all, which is
   a better outcome than restoring billing would have been.
2. **The study panel had to be regenerated**, because the 94 is pinned by
   `manifest/study_panel.tsv`, which lived in the dead bucket. `select_panel.py --every 2`
   reproduced it exactly: 95 datasets matching `rehearsal_binding_dinuc.csv`, of which 94 clear
   the GC floor and match `rehearsal_binding_gc.csv`.
3. **Three bugs, none in the new science.** A false `accelerator` field on non-CUDA devices; a
   test suite importing the sibling project's `rbp`; and a k-mer/SpliceBERT row-set mismatch
   caught by a composition-baseline equality check in the 7th decimal.
4. **The gate needed an arithmetic cross-check.** Every value assertion read
   `deep_contrast.csv` alone, so editing that one file passed all of them. Caught by attacking
   it, not by reading it.

---

## 7. The next experiment, and it is affordable: Horlacher's own negative sets

Checked 2026-08-31, because the genomics referee named this as the single highest-leverage
thing available and the whole question was whether the data exists.

**It does, and it is downloadable today.** But it is not cited in the paper: the data
availability statement names no repository, and the artifact is discoverable only through a
comment on GitHub issue #2.

| | |
|---|---|
| Zenodo | `10.5281/zenodo.10600977`, "Benchmark-RBP Samples", CC-BY-4.0 |
| File | `samples.tar.gz`, **379 MB**, md5 `c93a38b8bd684be2ccb5f6f82c6c4700` |
| Direct | `https://zenodo.org/api/records/10600977/files/samples.tar.gz/content` |
| Contents | **302 experiments x 5 CV folds**, each with `positive.bed`, `negative-1.bed`, `negative-2.bed` |
| Format | BED6 single-nucleotide crosslink sites, NOT sequences |
| Builds | **mixed**: ENCODE GRCh38, Mukherjee-PAR-CLIP and iONMF GRCh37 |
| Downloads to date | 33. Essentially undiscovered. |

**Why this is worth more than any further GPU spend.** Their negatives are sampled from
transcripts that contain a binding site, so they are **expression-controlled by construction**
-- which is the axis on which our own negative set is weaker than the prior art we cite. Running
our nested decomposition on their negative-1 and negative-2 sets would at one stroke:

  * remove the expression confound (the referee measured 41.8% of our GC-arm negatives sit in
    genes with TPM < 1, and a single scalar -- log TPM -- separates our classes at AUROC 0.833);
  * remove the dependence on our own sampler's undocumented `pool_multiple` / greedy assignment;
  * convert "we re-derived a known effect on our own windows" into "we show what the known
    effect does to measured model value, on the field's own published benchmark", at 302
    experiments rather than 94.

**Cost: $0 and roughly a day.** The paper's primary model class is a 4-mer logistic and the
composition baseline is 19 features -- both CPU. No GPU is needed to answer the question for
the model the headline is about. Extending it to the CNN would be free locally; SpliceBERT
would not fit the remaining budget and is not required.

**What the work actually is.** `bedtools getfasta` against two genome builds (GRCh38 is already
on this disk at `rna-binding-proteins/data/raw/`; GRCh37 would need downloading, or restrict to
the ENCODE subset and avoid the issue entirely), choose a window size to match ours (101 nt),
then run the existing `gain_over_composition` unchanged.

**Do NOT try to regenerate their pipeline.** Their raw ENCODE input is a symlink to a private
Helmholtz Munich HPC path, no accession list is published, the required GENCODE transcript BED
is gitignored, the top-level Snakefile references directories that do not exist, and the fold
split uses `sort --random-sort` with no seed -- so their CV folds are not reproducible from
scratch. The Zenodo tarball is the only way to match their splits, which is precisely why it is
the thing to use.

**Caveat to check before trusting it.** The release covers 302 experiments; the paper reports
313. The 11-dataset shortfall is unexplained and should be resolved or stated.
