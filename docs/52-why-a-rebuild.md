# 52. Why a rebuild, when the study already worked

The science was finished. Four results, verified, two adversarial attacks tested and
survived. So why spend a day rebuilding the machinery that produced them?

Because **"it worked" and "it is reproducible" are different claims**, and only the first one
was true.

---

## Part 1. The four things that made it unreproducible

Not opinions. Each is a specific, checkable defect.

### 1.1 Eighteen files hardcoded one GCP project

```bash
$ grep -rl "rbp-composition-2026" --include="*.py" --include="*.sh" .
18 files
```

This is invisible for as long as you only ever run in that project. It is a *total* failure
the first time anyone else runs it — and "anyone else" includes you, on a new project, in six
months, when a reviewer asks.

**Why it is worse than it looks.** A hardcoded project id does not fail loudly. It fails by
*writing into the wrong account* or *reading a stale bucket and finding nothing*. Consider:

```python
# the old code
bucket = storage.Client(project="rbp-composition-2026").bucket("rbp-composition-2026-derived")
```

Run that from a new project with the old credentials still cached and you have silently
appended to the previous study. Run it without access and you get a 403 that looks like a
transient permission blip. Neither says "this pipeline is not portable".

**The fix**, `src/rbp/utils/cloud.py`:

```python
def _resolve(explicit, env_key, cfg_key, what):
    if explicit:
        return explicit
    if os.environ.get(env_key):
        return os.environ[env_key]
    v = _from_config().get(cfg_key)
    if v:
        return v
    raise RuntimeError(f"{what} is not configured. ...")
```

Four things to notice, because each is a decision:

1. **Resolution order is most-specific-first.** An explicit argument beats the environment,
   which beats the config file. That is the order of increasing generality, so a caller can
   always override.
2. **There is NO default.** This is the important one. A default is how you write into
   someone else's bucket. If the environment is not configured, that is a bug in the
   environment and the program should stop.
3. **Buckets derive from the project id.** `f"{project()}-derived"`. Bucket names are
   globally unique across all of Google Cloud — not per project, *global*. So a project id
   nobody has used gives you bucket names nobody has used, for free, and one variable
   configures the whole system.
4. **`@lru_cache` on the config read**, because this is called in loops and reading YAML per
   call is silly.

And a test, so the habit cannot come back:

```python
@pytest.mark.parametrize("path", sorted(_files(), key=str), ids=lambda p: str(p.name))
def test_no_hardcoded_project_id(path):
    hits = [...]
    assert not hits, f"{path} hardcodes the project id: ..."
```

Parametrised per file rather than one test over all files, so a failure names the offender
instead of dumping a list.

### 1.2 The dataset panel was a command-line flag somebody typed once

The study runs on 95 datasets. Where did 95 come from? From this, typed during one sweep:

```bash
python scripts/cloud_train.py manifest --arm dinuc --every 2
```

`--every 2` on a 189-dataset panel keeps every second one: 95. Nothing recorded that choice.
The consequence was that the project carried **four different dataset counts** — 189, 187,
95, 94 — each correct for a different question and none written down, so every conversation
about them started from scratch.

**This is the single largest source of confusion in the whole study**, and it is an
architecture problem, not a communication problem. A number that exists only in shell history
is not a parameter, it is an accident.

**The fix** is to make the panel an *artefact*. `scripts/select_panel.py` computes it once,
uploads it to `manifest/study_panel.tsv`, and every later stage reads it:

```python
if blob.exists() and not a.force:
    log(f"panel already exists: {len(d)} datasets. Refusing to redefine it.")
    log("Pass --force only if you intend every downstream result to be invalidated.")
    return
```

It refuses to redefine itself, because silently redefining the panel mid-study makes half the
results describe a different population from the other half — and nothing anywhere would say
so.

It also **asserts its own validity**:

```python
lo_ok = picked.pairs.min() <= full.pairs.quantile(0.05)
hi_ok = picked.pairs.max() >= full.pairs.quantile(0.95)
if not (lo_ok and hi_ok):
    raise SystemExit("panel is size-biased: ... Refusing to write it.")
```

Why this specific check? Because AUROC correlates with dataset size at **r = +0.53 to +0.67**
across every model class. So "keep the biggest N" would confound the panel with the very
quantity being measured — the subset would look better than the population for a reason with
nothing to do with the biology. Systematic sampling by pair rank is unbiased in size *by
construction*, and the assertion proves the sample actually spans both tails.

### 1.3 Four stages ran only on a laptop

- variant assignment (needs the 3.1 GB genome)
- phyloP conservation fetch (needs UCSC over HTTP)
- aggregation of every result table
- figure rendering

Two of the paper's four results depended on files that existed on exactly one machine. The
final assembly — merging arms, building the four-model table, computing the cluster-corrected
ladder — was a sequence of ad-hoc interactive commands.

**That is the part of any study most likely to be irreproducible**, because it is the part
nobody writes down. It is done once, at the end, in a terminal, and the numbers in the paper
come from whatever was in memory that afternoon.

### 1.4 Nothing checked that a rerun produced the same science

There was no verification step. A pipeline that runs to completion and quietly produces
different numbers is **worse** than one that crashes, because nobody diffs a plausible table.

---

## Part 2. The architecture

### 2.1 The shape of it

```
                        YOUR LAPTOP
                   (submits and reads; never computes)
                              |
        +---------------------+---------------------+
        |                                           |
    GCP BATCH                                    MODAL
  (CPU fan-out)                              (GPU, no quota gate)
        |                                           |
   +----+----+----+----+                    +-------+-------+
   |    |    |    |    |                    |       |       |
 ingest panel prep rehearsal analysis   splicebert locality clinvar
   |    |    |    |    |                    |       |       |
   +----+----+----+----+--------------------+-------+-------+
                              |
                    GCS  (the only shared state)
        raw/            immutable inputs
        processed/      matched datasets, both arms
        panel/          per-arm candidate lists
        manifest/       task lists AND study_panel.tsv
        runs/           trained weights, per-fold scores
        variants/       ClinVar scores, three arms
        results/        tables and figures
```

**GCS is the only shared state, and that is the whole design.** No stage talks to another
stage. Every stage reads objects and writes objects. This is what makes the pipeline
resumable, parallelisable across two clouds, and debuggable — you can always answer "what did
this stage actually produce?" by listing a prefix.

### 2.2 Why two clouds

Not architectural elegance. A wall.

```
GPUS_ALL_REGIONS: limit 0     <- and the increase request is auto-denied,
                                 NOT_ENOUGH_USAGE_HISTORY, for 8, for 4, for 1
CPUS_ALL_REGIONS: limit 12
```

AWS returns 0 on all four GPU families. Azure forbids GPU quota on a free trial outright.
Modal gates nothing.

Measured, not assumed: transformer inference on the e2 CPUs GCP *will* give you is **4.9×
slower** than the laptop. Three nodes at 1/4.9 speed is slower than one Mac, so "run it on
GCP CPU" was never a real option for the transformer stages.

So: **CPU fan-out on Batch, GPU on Modal.** The split is forced by quota, and it turns out to
be free — a training task's only cloud dependency is `google-cloud-storage`, so the same
script runs unchanged on both platforms.

### 2.3 Why the stages are ordered the way they are

The order is not arbitrary and one edge is genuinely surprising.

```
1 terraform ──> 2 images ──> 3 ingest ──> 4 panel ──> 5 prep ──> 6 select
                                                                    |
                    +-----------------+-----------------+-----------+
                    |                 |                 |
                7 rehearsal      8 cnn (GPU img)   11 variants
                    |                 |                 |
                   R1                R2                 |
                                      |                 |
                              9 splicebert (Modal)      |
                                      |                 |
                            +---------+---------+       |
                            |                   |       |
                      10 locality          12 clinvar <-+
                            |                   |
                           R3                  R4
                            |                   |
                            +-------> 13 analysis ──> 14 verify
```

**The surprising edge: selection comes AFTER preprocessing.** I got this wrong in the first
plan and it would have failed on the fresh project immediately.

`pairs` — the number the panel is size-ranked on — counts the positives that could *actually
be matched* to a negative. Matching is a search. It succeeds to different degrees per
dataset. So `pairs` is a **result of preprocessing, not an input to it**. You cannot select a
size-ranked panel before prep has produced the counts.

Nothing is saved by trying to invert this: full prep over all 244 candidates × 2 arms is
about $2. The savings from a smaller panel come from the GPU stages, which are 95% of the
cost.

**The three-way parallel fan-out** after stage 7 is real and worth taking. Stages 8, 9 and 11
share no inputs and consume **different quota pools** — GCP vCPU versus Modal containers — so
they contend for nothing. Running them sequentially wastes the duration of the shorter two.
`./run.sh parallel` does this.

Stages 10 and 12 are genuinely downstream: locality needs SpliceBERT's fold-0 weights,
ClinVar needs those *plus* the variant tables.

---

## Part 3. Design decisions, with the alternative that was rejected

Every one of these is a place where a different choice was available and defensible.

| decision | alternative | why this one |
|---|---|---|
| **GCS as the only shared state** | a database, or passing files between stages | listing a prefix answers "what did this produce?" with no tooling. A database adds a service to run and a schema to migrate for data that is written once and read many times |
| **Completion marker written LAST** | write a marker first, then the payload | a task preempted between the two must **redo** its work, not be skipped. Marker-first means a crash produces a permanent silent gap |
| **No default project id** | fall back to a sensible one | a default is how you write into another account. Failing loudly costs one error message; guessing costs a corrupted study |
| **Task counts read from the manifest** | pass `--count N` | a count you typed is a count that will be wrong the first time the panel changes size. This exact bug failed a job whose every real task had succeeded |
| **Image pinned by digest** | pin by tag | a tag is mutable. "Which image produced this result?" is only answerable by digest, which is also why the model weights are baked into the image rather than downloaded at run time |
| **Weights baked into the image** | fetch from HuggingFace at run time | workers have Private Google Access and no route to `huggingface.co`. Baking makes weights and code share one digest. Cloud NAT was the alternative and costs per hour to give workers internet they should not have |
| **One generic submitter** | one script per stage | seven near-copies drift. The hardcoded-count bug existed in one copy and not the others |
| **Verification as a pipeline stage** | check by hand at the end | a hand check is not run when you are tired, and that is exactly when it matters |
| **Preflight refuses to spend** | let stages fail and read the error | every failure of the first build was a free-to-detect environment problem discovered *after* spending |
| **On-demand VMs by default** | spot everywhere | at 12 vCPU there is no headroom to replace a preempted worker. 7 vCPU-hours costs $0.07 on spot and $0.23 on demand — sixteen cents to remove a class of failure |
| **Parallelism 8, not 12** | use the full quota | VM creation fails **at** the limit, not approaching it. Asking for exactly the quota gets you an infinite retry loop and an alarming error that changes nothing |
| **Per-job service accounts** | one identity for all jobs | four accounts exist so a preprocessing task cannot write a model and an ingest task cannot touch results. One identity discards that for no gain |
| **Stage in / run existing code / stage out** | reimplement the logic against GCS | strand handling, window offsets and ref-allele checks are subtle. A second copy diverges silently from the copy that produced the published numbers. Copying bytes is cheap; copying logic is how two versions of a result appear |

---

## Part 4. What this bought

Running the rebuilt pipeline on a genuinely new project surfaced **sixteen bugs**, every one
of them invisible in the working project. Full accounting in
[55-the-bug-chronicle.md](55-the-bug-chronicle.md). The headline:

- one would have **deleted the original study's results** with no prompt
- one would have written a **zero-dataset panel** that every downstream stage trusted
- one would have failed **475 training tasks** on `import torch`
- one shipped a container **missing a source package**, caught only by the image's own tests

And three independent confirmations that the reproduction is real, not merely running: the
raw bucket totals **3.97 GiB** against the original's 3.97 GiB; the panel is **139 + 105**
against 139 + 105; and per-dataset pair counts match exactly (`SAFB2:K562 → 3,744 pairs` in
both).

That last one is the point of the whole exercise. Nothing else would have told us.
