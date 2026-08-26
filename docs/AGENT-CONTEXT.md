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

| | result |
|---|---|
| **R1** | GC-matched 0.796 → dinuc-matched 0.689, cost **−0.1070**, **187/187 datasets fall**, p=1.9e-32. **The finding is not the drop**: gain over composition **+0.027 → +0.064**, it DOUBLES |
| **R2** | composition 0.628 < k-mer 0.688 < CNN 0.708 < SpliceBERT 0.809, identical splits, 95 datasets / 79 proteins |
| **R3** | ISM Gini: SpliceBERT more positionally concentrated on **89/95**, reversed on **0/95**, median +0.062, p=3.9e-17 |
| **R4** | ClinVar ladder, cluster-corrected: k-mer 0.529 → **wrong-protein head 0.656** → right-protein head 0.829; conservation 0.906. Coefficients 0.102 / 0.622 / **1.616** over 1,426 genomic blocks |

**R4's honest reading.** A wrong-protein head already beats the k-mer, so part of the signal
is generic sequence plausibility inherited from pretraining. Binding-specific signal is the
GAP between rungs, not the top number. Conservation still beats all of it.

**The four panels** (`docs/PANELS.md` in both repos is the single source of truth):
`VARIANT 94 ⊂ DEEP 95 ⊂ FULL 189`, and `MATCHED 187 ⊂ FULL 189`. Never mix them without
saying so.

## 3b. RESULTS REPRODUCED on rbp-repro-2026 (2026-08-25)

Three of four, all inside golden.yaml tolerances (8/8 checks pass).

| | reproduced | reference |
|---|---|---|
| **R1** | -0.1095, **94/94 fall**, p=3.81e-17, gain +0.0154 -> +0.0607 (**3.94x**) | -0.107, 3.61x |
| **R2** | composition 0.6279, k-mer 0.6875, CNN 0.7063, SpliceBERT 0.8091 | 0.628/0.688/0.708/0.809 |
| **R3** | median +0.064, more local 91/95, **significantly reversed 0/95**, p=3.94e-17 | +0.062, 91/95 |
| **R4** | OUTSTANDING -- blocked on stage 11 | |

Integrity: 0 NaN, 0 at chance, no duplicates, all 95 datasets have full 5-fold sets for both
CNN and SpliceBERT. Sweeps: rehearsal 189/189, CNN 475/475, SpliceBERT 475/475, locality 95/95.

**R1's n is 94, not 187** -- restricting the study to 95 datasets means only 94 clear the
min_pairs=400 floor in BOTH arms (NCBP2:K562 matches 384 pairs under GC, 406 under dinuc).
Stage 7 is 189 TASKS: 95 dinuc + 94 gc, one per dataset per arm.

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
