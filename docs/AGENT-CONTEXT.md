# AGENT CONTEXT — read this first if you are an assistant resuming this work

Compressed, load-bearing state. Written 2026-08-25. If the conversation was compacted, this
file plus `git log` is enough to continue without re-deriving anything.

---

## 1. Two repositories, and they are not the same thing

| path | what it is | status |
|---|---|---|
| `Deep Learning Project/rna-binding-proteins/` | **the original study.** All four results live here, verified. Do not break it. | scientifically COMPLETE |
| `Deep Learning Project/rbp-repro/` | **the reproducible rebuild.** Same science, rebuilt as cloud stages. | mid-run, stage 5 |

The original's results are the ground truth the rebuild is checked against
(`rbp-repro/config/golden.yaml`). **A near-miss deleted them once** — see §6.

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

## 4. Where the rebuild is right now

```
0 preflight   DONE   21 passed
1 terraform   DONE   78 resources, 4 buckets, 6 service accounts
2 images      DONE   cpu sha256:9265f8c8..., gpu sha256:23e7d9cb...
3 ingest      DONE   3.97 GiB, byte-identical total to the original
4 panel       DONE   139 K562 + 105 HepG2, identical to the original
5 prep        RUNNING  488 tasks, ~2.2/min, ETA ~16:00, job prep-0825-121949
6 select      next   writes manifest/study_panel.tsv -- THE panel
7 rehearsal   -> R1
8 cnn         GPU image (it has torch; the cpu image does not)
9 splicebert  MODAL, ~$31, ASK THE USER FIRST -- this is a standing instruction
10 locality   MODAL
11 variants   needs a public IP (UCSC phyloP)
12 clinvar    MODAL
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
