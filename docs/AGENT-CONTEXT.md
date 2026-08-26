# AGENT CONTEXT — read this first if you are an assistant resuming this work

Compressed, load-bearing state. Written 2026-08-25. If the conversation was compacted, this
file plus `git log` is enough to continue without re-deriving anything.

---

## 1. Two repositories, and they are not the same thing

| path | what it is | status |
|---|---|---|
| `Deep Learning Project/rna-binding-proteins/` | **the original study.** All four results live here, verified. Do not break it. | scientifically COMPLETE |
| `Deep Learning Project/rbp-repro/` | **THE PROJECT.** The first build is being discarded; only this one is published or discussed. | stages 0-10 done, 11 running |

**The user has decided the first build is discarded entirely** -- not published, not
mentioned, treated as a standalone project in interviews. So there is no "two swapped
proteins" story and no "reproduce the original" framing: `golden.yaml` holds THIS
pipeline's reference values. The first build survives only as a cross-check while the rebuild
finishes. **A near-miss deleted it once** -- see §6.

## 2. Two GCP projects, two Modal accounts

| | project | notes |
|---|---|---|
| original | `rbp-composition-2026` | holds the 485 SpliceBERT checkpoints, ~$31 of GPU time |
| rebuild | `rbp-repro-2026` | created 2026-08-25, same billing account `017994-4FC1D0-8D5176` |
| Modal (old) | default profile | $30 credit exhausted, ~$6 out of pocket spent |
| Modal (new) | profile `meghanasai1802` | $30 credit, secret `rbp-gcp` created and verified |

Buckets are always `{project}-derived`, `{project}-raw`, `{project}-artifacts`,
`{project}-tfstate`. Nothing hardcodes a project: everything resolves through
`src/rbp/utils/cloud.py`, and `tests/unit/test_no_hardcoded_project.py` fails the build if a
literal reappears.

## 3. The science, in numbers you must not misquote

**Central claim.** Negative-set construction, not architecture, determines what an RBP
binding benchmark measures — and under proper matching the model earns *more* credit, not
less.

| | result, as of 2026-08-26 |
|---|---|
| **R1** | Composition share: a 19-feature composition model recovers **94.8%** [92.1, 97.4] of the k-mer model's skill above chance under GC matching, **67.8%** [62.1, 73.8] under dinucleotide matching. Drop **27.0 points** [23.1, 30.9]. Nested gain over composition +0.0265 -> +0.0662 (**2.50x**), significant on 85-87% of datasets. AUROC 0.798 -> 0.689, cost -0.1095, **94/94 fall**, p=3.81e-17 |
| **R2** | composition 0.6279 < k-mer 0.6875 < CNN 0.7063 < SpliceBERT 0.8091, identical splits, 95 datasets / 79 proteins. **Not a finding** -- Horlacher 2023 covers this with 11 methods. Methods table only |
| **R3** | ISM Gini: SpliceBERT more concentrated on **91/95**, significantly reversed on **0/95**, median +0.064, p=3.94e-17. **Cut from the paper**: linear models give diffuse attributions and transformers spiky ones by construction, so it measures model class, not biology |
| **R4** | **PAIRED per dataset, not pooled.** At >=20 pathogenic variants (n=44): right-protein head **0.7553** vs wrong-protein **0.6908**, gap **+0.0645**, 33/44, p=3.9e-04. At >=50 (n=31): +0.1031, 27/31. Conservation **0.8921** leads everything. Within-dataset conservation-controlled coefficients **0.716** [0.640, 0.801] vs **0.445** [0.365, 0.523], non-overlapping |

### THREE CORRECTIONS MADE ON 2026-08-26. Do not reintroduce the old numbers.

**1. R4's pooled ladder was inflated and is no longer the reported result.** Pooling ~19k
variants into one AUROC per arm partly measures WHICH DATASET a variant came from: mean
|delta| per dataset correlates with that dataset's pathogenic rate at Spearman **+0.73** and
spans **10.4x**. Matched pooled 0.829 against 0.755 paired; the gap +0.149 is really +0.065.
Conservation was the only arm immune, because phyloP is on a fixed external scale -- which is
why it stayed invisible: the arm that could not be inflated was winning anyway.

**2. A TRIVIAL POSITIONAL BASELINE BEATS THE MODEL.** "What fraction of the OTHER variants in
this 1-Mb window are pathogenic", leave-one-out, no sequence and no model, reaches **0.8139**
where the model reaches 0.7553. The model does not beat it in any stratum: 15/44 wins at >=20
pathogenic (**p=0.007 against us**), a tie at >=50 (p=0.53). **So absolute AUROC on
peak-proximal ClinVar variants is uninformative about model utility.** No version of this
work may report 0.755 as evidence of usefulness. What survives is the PAIRED specificity
contrast, which the baseline cannot explain because it applies equally to both arms.

**3. R1's gain was being verified as the wrong quantity.** The gate computed
(auroc - composition_auroc), a difference of two separately fitted models with no interval
and no p-value, giving 3.94x. The claim is the NESTED gain -- composition alone versus
composition plus the score -- which the rehearsal already reports as delta_auroc with a CI
and a p-value: **2.50x**. Smaller and defensible.

**The wrong-protein control was attacked and held.** Three checks, all negative for
contamination: the floor does not depend on the donor sharing a cell line (p=0.83), only
weakly tracks the donor's own strength (rho=+0.23), and stays flat at ~0.69 across power
strata while the matched arm climbs 0.66 -> 0.80. It is a generic-plausibility floor, not a
similarity artefact.

**The panels.** `docs/PANELS.md` is the single source of truth. This pipeline: study panel
**95**, GC arm **94** (NCBP2:K562 matches 384 pairs under GC and 406 under dinucleotide, so
it clears the 400 floor in one arm only), stage 7 is **189 TASKS** = 95 + 94, variant arm
**95** with **82** usable for per-dataset AUROC and **44** adequately powered. The earlier
study's 187/DDX51 framing belongs to the discarded build and must not reappear.

## 3b. VERIFICATION STATUS

**56/56 golden checks pass** reading from GCS; 575 unit tests pass. Integrity: 0 NaN, 0 at
chance, no duplicates, all 95 datasets have full 5-fold sets for both CNN and SpliceBERT.
Sweeps: rehearsal 189/189, CNN 475/475, SpliceBERT 475/475, locality 95/95, ClinVar 95/95
matched and 95/95 mismatched.

**The gate itself failed twice today and both are recorded in golden.yaml.** It passed 33/33
while certifying an inflated pooled statistic, and separately certified a 3.94x gain ratio
that nothing in the write-up claimed. Same failure both times: checking the number the code
produced rather than whether that number was the right number. The unflattering results --
the all-datasets stratum that shows nothing, and `model_minus_prevalence_max` -- are now
asserted so they cannot be quietly dropped.

## 4. Where the rebuild is right now

```
0 preflight   DONE   21 passed
1 terraform   DONE   78 resources, 4 buckets, 6 service accounts
2 images      DONE   REBUILD AFTER EVERY CODE CHANGE (see 5b)
3 ingest      DONE   3.97 GiB, byte-identical total to the original
4 panel       DONE   139 K562 + 105 HepG2, identical to the original
5 prep        DONE   488/488, panel counts identical to reference
6 select      DONE   95 datasets, 79 proteins, deterministic sort
7 rehearsal   DONE   189/189 -> R1 REPRODUCED
8 cnn         DONE   475/475 -> R2 REPRODUCED (gpu image + --device cpu)
9 splicebert  DONE   475/475 on FREE credit, $0 out of pocket
10 locality   DONE   95/95 -> R3 REPRODUCED
11 variants   RUNNING after 3 fixes (SA, raw read, get_blob, all peaks)
12 clinvar    NEXT, ~$0.60 Modal -> R4
13 analysis
14 verify     asserts golden.yaml; exit 1 means the science did not reproduce
```

Run with: `export GOOGLE_CLOUD_PROJECT=rbp-repro-2026; ./run.sh stage N`.
`PY=../rna-binding-proteins/.venv/bin/python` and
`PATH=../rna-binding-proteins/.venv/bin:$PATH` (modal must be on PATH, not just installed).

## 5. Standing user instructions — these persist

- **Ask before every spend**, and specifically before anything on Modal.
- Modal budget: **$20 out of pocket + $30 credit**. Stage 9 is ~$31 and is 95% of it.
- No local compute in the rebuild. Everything on Batch or Modal.
- Documentation must teach cloud **three ways**: CLI/Terraform, UI step by step, and backend
  internals. Every file, every line, every bug.
- No em-dashes in the user's own deliverables (paper, slides). Code stays concise, minimal
  comments, no echo spam.
- The user is testing the pipeline on a fresh project deliberately. **Bugs found are the
  point, not a setback.** Report them plainly.

## 5b. THE FOUR TRAPS THAT COST THE MOST TIME (full record: docs/58-the-run-chronicle.md)

**A stale image makes a job succeed while doing nothing.** The rehearsal reported 189/189
SUCCEEDED with an empty rehearsal/gc/, because the container predated the arm-from-row change
and every task took ARM=dinuc from the environment, hit an existing marker, and exited 0.
After ANY code change: `gsutil cat gs://${PROJECT}-artifacts/images/cpu_digest.txt`.

**The laptop keeps sneaking back into the pipeline, on both clouds.** Once as a shell loop
sequencing two Batch jobs; three hours later as `modal run` without `--detach`, which ties the
app's lifetime to the local client and killed 475 GPU tasks at 43%. Fixed the pattern on GCP
and did not think to look for it on Modal.

**A stage-in list is a contract with every function the stage calls.** Stage 11 failed three
times: wrong service account, then missing read on raw, then missing ENCODE peaks -- twice,
because I first staged peaks only for the study panel and `--what assign` walks the full
candidate panel. 19.7 MiB. Being selective saved nothing.

**Safety mechanisms fire on wrong instructions, and that is them working.** An IAM condition
(twice), a device guard, a completion marker, an app-lifetime default. Five of the 24 bugs.
Each cost minutes instead of a corrupted result. Do not weaken the guard; fix the instruction.

## 6. The mistakes worth remembering, because they will recur

**The one that nearly cost everything.** `main.tf` hardcoded the Terraform backend bucket to
`rbp-composition-2026-tfstate`. Pointing the config at the new project made `terraform plan`
load the OLD project's state and propose `63 to add, 1 to change, 63 to destroy` — the 63
destroys being the original study's buckets. `run.sh` ran `apply -auto-approve`. Caught by
reading the plan. Now: partial backend config passed at init with `-reconfigure`, plus a
guard that counts destroy actions and refuses to apply if any exist.

**Diagnosing from the story you already have.** I called a vCPU drop "spot preemption",
switched to on-demand, and measured no improvement. The job's events said
`CODE_GCE_QUOTA_EXCEEDED` and contained zero preemption events. The real cause: parallelism
12 × 4 tasks/node asks for exactly 12 vCPU, which is exactly the quota, and **VM creation
fails at the limit rather than reaching it**. Read the evidence before choosing the cause.

**Guards that hide other failures.** `figures.py` skips missing tables quietly so it can run
mid-pipeline. That is correct — and it silently hid two figures reading table names nothing
wrote. A tolerant gate needs a strict counterpart somewhere.

**Tests only certify where they look.** `test_no_hardcoded_project.py` searched `scripts`,
`src`, `cloud` — not `docker/`. Two hardcoded project ids survived in the cloudbuild files,
in a test written to prevent exactly that.

**Estimating from one sample.** I quoted 45–50 minutes for prep from a single 47.5s task,
early in a biggest-first manifest. Measured over a window: 2.2/min, 3.5 hours. Measure.

## 7. What is NOT done

- The **manuscript**. Not one word. Every number exists in committed tables; realistically
  2–3 days of writing.
- Stages 6–14 of the rebuild.
- `docs/53-cloud-from-zero.md` and the rest of this documentation set are in progress.
- Modal token/secret pasted in chat on 2026-08-25 — **remind the user to rotate it** when the
  run finishes.
