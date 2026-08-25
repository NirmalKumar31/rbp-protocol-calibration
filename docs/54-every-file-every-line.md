# 54. Every file, every line

The complete code walkthrough. For each file: what it is for, how it is invoked, what every
non-obvious line does, and why it is written that way rather than the obvious alternative.

Files are grouped by role, not alphabetically, because the roles are what you need to hold in
your head.

```
rbp-repro/
├── run.sh                      THE ORCHESTRATOR. 15 stages.
├── cloud/
│   ├── submit.sh               one Batch submitter for every stage
│   ├── terraform/*.tf          all infrastructure
│   └── modal/                  the two Modal apps
├── scripts/
│   ├── preflight.py            stage 0: refuse to spend on a broken environment
│   ├── select_panel.py         stage 6: define THE panel, once
│   ├── cloud_ingest.py         stage 3: raw inputs
│   ├── cloud_prep.py           stage 5: matched datasets
│   ├── cloud_rehearsal.py      stage 7: composition + k-mer  -> R1
│   ├── cloud_train.py          stages 8/9: CNN and SpliceBERT -> R2
│   ├── locality_probe.py       stage 10: ISM Gini             -> R3
│   ├── cloud_variants.py       stage 11: ClinVar + phyloP
│   ├── variant_splicebert.py   stage 12: the ladder           -> R4
│   ├── cloud_analysis.py       stage 13: tables + figures
│   ├── figures.py              the five figures
│   └── verify.py               stage 14: assert golden.yaml
├── src/rbp/                    the library. Platform-agnostic.
├── config/
│   ├── params.yaml             every parameter, including cloud:
│   └── golden.yaml             what a correct reproduction must produce
└── tests/                      480 tests
```

---

# Part 1. `run.sh` — the orchestrator

## 1.1 Header and resolution

```bash
set -uo pipefail
```

`-u` errors on an unset variable, which catches typos in variable names. `-o pipefail` makes a
pipeline fail if any stage fails, not just the last.

**`-e` is deliberately absent.** With `-e`, a failing stage kills the whole script. That is
right for CI and wrong here: the point of an overnight run is to come back to as much finished
work as possible **plus an honest log of what broke**, not a script that died at step 2.

```bash
PROJECT_ID="${GOOGLE_CLOUD_PROJECT:-$($PY -c '...cloud.project()')}"
DERIVED="${DERIVED_BUCKET:-${PROJECT_ID}-derived}"
export GOOGLE_CLOUD_PROJECT="$PROJECT_ID" DERIVED_BUCKET="$DERIVED" RAW_BUCKET="$RAW"
```

The environment wins over the config file; the config file is the fallback. The `export` is
what makes every child process — python, gsutil, submit.sh — see the same answer. Without it,
`run.sh` and the scripts it calls could disagree about which project they are in, which is
exactly the class of bug this whole design exists to prevent.

## 1.2 The paid gate

```bash
confirm() {
  local what=$1 cost=$2
  say "PAID STAGE: $what  (estimated $cost)"
  if [ "${RBP_YES:-0}" = "1" ]; then say "RBP_YES=1, proceeding"; return 0; fi
  printf 'Type YES to spend %s on "%s": ' "$cost" "$what"
  read -r ans; [ "$ans" = "YES" ] || die "not confirmed"
}
```

Three deliberate choices. It requires the literal string `YES`, not `y`, because `y` is muscle
memory. It **prints the estimate**, so the decision is informed. And `RBP_YES=1` exists so the
whole thing can run unattended *once you have decided* — a gate you cannot bypass gets worked
around with something worse.

## 1.3 Two gates, not one

```bash
gate_preflight() { [ -f .preflight-ok ] || die "run ./run.sh preflight first"; }
gate_modal()     { [ -f .preflight-modal-ok ] || die "..."; }
```

These are separate because a single all-or-nothing preflight is **unsatisfiable on a fresh
project**. The `rbp-gcp` secret holds a key for a service account that does not exist until
stage 1 has run. So a combined gate would demand a secret that cannot yet be created, in order
to permit the stage that creates its prerequisite. GCP stages need the GCP gate; Modal stages
need both.

## 1.4 The stages

```bash
STAGES=(s0_preflight s1_terraform s2_images s3_ingest s4_panel s5_prep s6_select \
        s7_rehearsal s8_cnn s9_splicebert s10_locality s11_variants s12_clinvar \
        s13_analysis s14_verify)
```

Index *n* maps to `s{n}_*`, so `./run.sh stage 7` runs `s7_rehearsal`. Verified in bash: 15
elements, all aligned. (My own audit script once reported this as broken because its regex
treated the `\` line continuations as array elements — the script was wrong, not `run.sh`.)

### `s1_terraform` — the destroy guard

```bash
local TFSTATE="${PROJECT_ID}-tfstate"
if ! gsutil ls -b "gs://${TFSTATE}" >/dev/null 2>&1; then
  gcloud storage buckets create "gs://${TFSTATE}" ... --uniform-bucket-level-access
  gcloud storage buckets update "gs://${TFSTATE}" --versioning || true
fi
```

The state bucket is created **outside Terraform**, because storing state about the bucket that
stores the state is a bootstrap paradox. `--versioning` means a corrupted state file can be
rolled back — cheap insurance on an object written rarely.

```bash
terraform init -input=false -reconfigure -backend-config="bucket=${TFSTATE}"
```

`-reconfigure` is load-bearing. Without it, init reuses whatever backend a previous checkout
cached in `.terraform/`, which is precisely how another project's state leaks in and produces a
63-destroy plan.

```bash
terraform plan -input=false -no-color -out=tfplan.new
DESTROYS=$( terraform show -no-color tfplan.new | grep -cE "^  # .* will be destroyed" || true )
if [ "${DESTROYS:-0}" -gt 0 ]; then die "plan contains ${DESTROYS} DESTROY actions..."; fi
terraform apply -input=false -auto-approve tfplan.new
```

Four things:

1. `-out=tfplan.new` saves the plan, and `apply tfplan.new` applies **that saved plan**. Without
   the file, `apply -auto-approve` re-plans and can do something different from what you
   inspected.
2. The destroy count is the guard. A first apply on an empty project is additive **by
   definition**, so any destroy proves the state describes a different project.
3. `|| true` on the grep, because `grep -c` exits 1 when it finds nothing, and with `pipefail`
   that would abort the script on the *good* path.
4. `${DESTROYS:-0}` because an empty string in `-gt` is a syntax error.

### `s2_images` — concrete substitutions

```bash
local REPO="${REGION}-docker.pkg.dev/${PROJECT_ID}/rbp"
local SHA; SHA=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)
for kind in cpu gpu; do
  gcloud builds submit --config="docker/cloudbuild.${kind}.yaml" \
    --substitutions="_IMAGE=${REPO}/${kind},_ARTIFACTS=${PROJECT_ID}-artifacts,_GIT_SHA=${SHA}" .
done
```

`local SHA; SHA=$(...)` on two statements, not `local SHA=$(...)`: `local` always returns 0, so
combining them hides a failure in the command substitution.

Substitutions are passed explicitly because Cloud Build does **not** expand `$PROJECT_ID`
inside a user-defined substitution's *default value* — the literal string reaches docker and is
rejected for containing capitals.

### `s5_prep` — the forced ordering

```bash
$PY scripts/cloud_prep.py index    || die "prep index"
$PY scripts/cloud_prep.py manifest || die "prep manifest"
./cloud/submit.sh prep             || die "prep"     # one job, both arms
for arm in dinuc gc; do $PY scripts/cloud_prep.py finalize --arm "$arm"; done
```

Prep runs on **every candidate**, before the panel exists, and that order is forced: `pairs`
counts the positives that could actually be matched to a negative, so it is a *result* of
preprocessing. A size-ranked panel cannot be selected before those counts exist.

One job for both arms, because each task reads its arm from its own manifest row
(`cloud_prep.py:173`). Passing `--arm` to the job would have been wrong.

The `index` and `manifest` steps were **missing entirely** in the first version — the submitter
would have read a manifest that did not exist.

### `s9_splicebert` — the expensive one

```bash
gate_preflight; gate_modal
say "This is 95% of the money. Confirm your Modal balance is >= \$35 before continuing."
confirm "SpliceBERT, 5 folds per dataset, Modal A10G" "~\$31 OUT OF POCKET"
```

The words OUT OF POCKET are in the prompt because GCP credit and real money feel identical in
a terminal and are not.

### The parallel track

```bash
parallel) gate_preflight; gate_modal
          ( s8_cnn        > logs/track-cnn.log      2>&1; say "track CNN done rc=$?" ) &
          ( s9_splicebert > logs/track-sb.log       2>&1; say "track SpliceBERT done rc=$?" ) &
          ( s11_variants  > logs/track-variants.log 2>&1; say "track variants done rc=$?" ) &
          wait ;;
```

Three subshells, each with its own log, then `wait`. Safe because the three stages share no
inputs and consume **different quota pools** — GCP vCPU versus Modal containers. Stages 10 and
12 are excluded: locality needs SpliceBERT's fold-0 weights, ClinVar needs those plus the
variant tables.

---

# Part 2. `cloud/submit.sh` — one submitter for every stage

Replaces seven near-identical scripts. Seven copies drift, and the hardcoded-count bug existed
in one copy and not the others.

## 2.1 Task counts, never typed

```bash
manifest_rows() {
  local n
  n=$(gcloud storage cat "gs://${DERIVED}/${key}" 2>/dev/null | tail -n +2 | wc -l | tr -d ' ')
  [ -n "$n" ] && [ "$n" -gt 0 ] || { echo "EMPTY_MANIFEST"; return 1; }
  echo "$n"
}
```

`tail -n +2` skips the header. `tr -d ' '` strips the leading whitespace `wc -l` emits on macOS,
which otherwise breaks the numeric comparison. The empty check catches a manifest that exists
but has no rows — a stage that ran but produced nothing.

## 2.2 The per-stage table

Each stage sets: the service account, the script and arguments, the task count, parallelism,
tasks per node, machine type, per-task CPU and memory, whether it needs an external IP, disk
size and timeout.

```bash
prep)
    SA=rbp-prep
    SCRIPT="scripts/cloud_prep.py"; ARGS="prep"
    COUNT=$(manifest_rows "manifest/prep_tasks.tsv") || exit 1
    PAR=8; PER_NODE=4; MACHINE=e2-standard-4; CPU=900; MEM=3500
    EXTERNAL=0; DISK=100; TIMEOUT=7200 ;;
```

**`PAR=8`, not 12.** `CPUS_ALL_REGIONS` is 12 and `e2-standard-4` is 4 vCPU, so parallelism 12
with 4 tasks per node asks for exactly three nodes and exactly 12 vCPU. VM creation then fails
with `CODE_GCE_QUOTA_EXCEEDED`: the limit is a ceiling you cannot touch, not one you can reach.
Batch retries the third node forever, the job runs 8-wide anyway, and the events fill with an
error that looks alarming and changes nothing.

**`prep_tasks.tsv`, not `study_panel.tsv`.** Prep runs before the panel exists. Pointing this at
the study panel made the fresh run unrunnable at stage 5.

**`CPU=900` on a 4000-milli machine** gives 4 tasks per node. It also caps each task below one
core, which is why `OMP_NUM_THREADS=1` matters — see §2.6.

## 2.3 Which image

```bash
IMAGE_KIND=cpu
[ "$JOB_TYPE" = "sweep" ] && IMAGE_KIND=gpu
```

Only the sweep needs torch, and the CPU image has none by design. The CNN is a torch model, so
a sweep task on the CPU image would fail on `import torch` 475 times. The GPU image runs fine
on a CPU machine; it just carries CUDA it will not touch, which is the cheaper mistake.

## 2.4 Pinned by digest

```bash
DIGEST=$(gcloud storage cat "gs://${PROJECT}-artifacts/images/${IMAGE_KIND}_digest.txt" | tr -d '[:space:]')
IMAGE="${REGION}-docker.pkg.dev/${PROJECT}/rbp/${IMAGE_KIND}@${DIGEST}"
```

A tag is mutable; a digest is content-addressed. "Which image produced this result?" is only
answerable by digest. `tr -d '[:space:]'` because a trailing newline inside an image reference
is a parse error with a baffling message.

## 2.5 The network block

```bash
if [ "$EXTERNAL" = "1" ]; then
  NETWORK=""
else
  NETWORK='"network": {"networkInterfaces": [{"network": "projects/'"${PROJECT}"'/global/networks/rbp-net", "subnetwork": "projects/'"${PROJECT}"'/regions/'"${REGION}"'/subnetworks/rbp-workers", "noExternalIpAddress": true}]},'
fi
```

Nesting **and** names both mattered. Batch wants `allocationPolicy.network.networkInterfaces`,
one level deeper than the obvious guess, and rejects the flat form outright. And the network is
`rbp-net`, not `default`: the whole point of `network.tf` is that workers sit on a subnet with
Private Google Access and no route to the internet. Naming `default` would have *worked*, which
is what made it the dangerous half of that bug.

The trailing comma is inside the string because the external case emits nothing, and a dangling
comma in the JSON template would be a syntax error.

The quoting `'"${PROJECT}"'` closes the single-quoted string, inserts the expanded variable,
and reopens it — the standard idiom for interpolating into single quotes.

## 2.6 Thread pinning

```json
"environment": {"variables": {
  "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
  "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1"
}}
```

numpy and scipy size their thread pools from the cores they can **see** — the host's count, not
the container's cgroup limit. Four tasks on a 4-core node at 0.9 CPU each means four processes
each spawning four threads: sixteen threads over four cores.

Measured: NCBP2, a 406-pair dataset, computes in **7 seconds** on a laptop. Unpinned in a
container it ran **over 20 minutes without finishing**. Pinned: **26.7 seconds**. A 189-task
sweep would have looked like a hang rather than an error. It also fixes summation order inside
BLAS, which makes the numerics reproducible run to run.

## 2.7 Waiting

```bash
while :; do
  STATE=$(gcloud batch jobs describe "$JOB" --format="value(status.state)" 2>/dev/null)
  case "$STATE" in
    SUCCEEDED) exit 0 ;;
    FAILED|CANCELLED) exit 2 ;;
    "") echo "cannot read job state; retrying" ;;
    *)  echo "[$(date '+%H:%M:%S')] ${STATE}  ${COUNTS}" ;;
  esac
  sleep 60
done
```

`gcloud batch jobs submit` returns on **acceptance of the JSON**. Without this loop, stage 5's
next line would finalize the panel against an empty bucket.

**Exit 2 for FAILED, not 1**, because a failed job is not necessarily a failed run: preemption
leaves resumable failures, and a task count past the end of a manifest once failed a job whose
every real task succeeded. The caller decides.

The empty-string case **assumes alive** and retries. An earlier monitor on the original project
treated an unparseable response as completion and reported a run finished at 273 of 475 tasks.

---

# Part 3. `scripts/preflight.py` — stage 0

The most valuable file here, because every failure of the first build was one of these and every
one was found *after* spending.

## 3.1 The check helper

```python
def check(name, ok, detail="", fatal=True):
    status = "PASS" if ok else ("FAIL" if fatal else "WARN")
    results.append((status, name, detail))
    print(f"  [{status}] {name}" + (f"  {detail}" if detail else ""), flush=True)
    return ok
```

Three levels, not two. Some findings are advisory (docker absent — Cloud Build is used anyway)
and some are fatal (no billing account). Collapsing them would either block on nothing or let a
real problem through. `flush=True` so output appears in order when redirected.

```python
def sh(cmd, timeout=90):
    try:
        p = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
        return p.returncode, (p.stdout or "") + (p.stderr or "")
    except Exception as e:
        return 1, str(e)
```

Never raises. A broken tool is a **finding**, not a crash — preflight's job is to report
everything wrong in one pass, not stop at the first.

## 3.2 The emptiness trap

```python
rc, out = sh(f"gcloud projects describe {proj} --format='value(projectId)'")
check("project exists and is visible", rc == 0 and proj in out, ...)
```

`gcloud <noun> list` with no project returns an **empty list, not an error**, so a missing scope
reads exactly like a missing result. `describe` on a named resource fails loudly. This is why the
check is written against `describe` rather than `list`.

## 3.3 The budget check that matters

```python
bad = [b.get("displayName", "?") for b in budgets
       if (b.get("budgetFilter", {}) or {}).get("creditTypesTreatment") == "INCLUDE_ALL_CREDITS"]
check("budgets exclude credits", not bad, ..., fatal=not not bad)
```

`INCLUDE_ALL_CREDITS` subtracts free credit from reported spend, so on a project with $300 of
trial credit **reported spend is $0** until the credit is gone and every threshold is
unreachable. The budget exists, the alerts are configured, and nothing can ever fire.

`(b.get("budgetFilter", {}) or {})` handles both a missing key and an explicit `null`, which the
API returns for budgets with no filter.

## 3.4 Reporting GPU quota without failing on it

```python
gpus = q.get("GPUS_ALL_REGIONS", {}).get("limit", 0)
check("GPUS_ALL_REGIONS (informational)", True,
      f"limit={gpus:g}" + (" -- expected on a new project; stage 8 uses Modal" if gpus == 0 else ""),
      fatal=False)
```

Always passes, always prints. GPU quota **will** be 0 on a new project and cannot be raised, so
failing on it would block a pipeline that is designed around that fact. Reporting it stops
somebody spending a day trying to fix it.

Note it is **not filtered to `limit > 0`**. An earlier quota survey on the original project did
filter that way and thereby hid `GPUS_ALL_REGIONS = 0` from its own output — twice, leading me
to tell the user GPUs were available.

## 3.5 The bucket check that distinguishes three cases

```python
if rc == 0:
    check(f"gs://{b}", True, "exists and is readable by us")
elif "does not have storage.buckets.get" in out or "403" in out:
    check(f"gs://{b}", False, "name is taken by another project. Bucket names are global...")
else:
    check(f"gs://{b}", True, "does not exist yet; terraform will create it", fatal=False)
```

Bucket names are globally unique, so a **403 means somebody else owns this name** — a fatal,
unfixable-by-retry condition. A 404 means it does not exist yet, which on a fresh project is
expected and fine. Conflating them would either block a normal first run or let Terraform fail
later with a confusing message.

## 3.6 The human gate

```python
check("Modal credit confirmed by a human", False,
      "CANNOT be checked automatically. Open modal.com and confirm >= $35 available "
      "BEFORE stage 8. Re-run with --modal-credit-ok once verified.", fatal=True)
```

Modal exposes no balance via CLI, and stage 9 is 95% of the spend. Rather than assume, this
fails until a human says otherwise. An honest unfixable check beats a silent assumption.

---

# Part 4. `scripts/select_panel.py` — stage 6

## 4.1 Refusing to redefine

```python
if blob.exists() and not a.force:
    log(f"panel already exists: {len(d)} datasets. Refusing to redefine it.")
    log("Pass --force only if you intend every downstream result to be invalidated.")
    return
```

Silently redefining the panel mid-study makes half the results describe a different population
from the other half, with nothing anywhere saying so. The message names the consequence rather
than just refusing.

## 4.2 Systematic sampling, and proving it

```python
picked = full.iloc[::every].reset_index(drop=True) if every > 1 else full
```

`full` is sorted by `pairs`, so `[::2]` keeps every second dataset **by size rank**. Unbiased in
size by construction: the sample spans the full range with the same shape.

The alternative — keeping the largest N — is confounded, because AUROC correlates with dataset
size at r = +0.53 to +0.67. A size-selected subset would look better than the population for a
reason with nothing to do with biology.

```python
lo_ok = picked.pairs.min() <= full.pairs.quantile(0.05)
hi_ok = picked.pairs.max() >= full.pairs.quantile(0.95)
if not (lo_ok and hi_ok):
    raise SystemExit("panel is size-biased: ... Refusing to write it.")
```

The assertion proves the claim rather than asserting it. If systematic sampling is ever broken —
by a sort that silently fails, say — the panel is not written.

## 4.3 Reporting the other arm

```python
both = set(picked.dataset) & set(oth.dataset)
log(f"of these, {len(both)} also clear the floor in the {other} arm "
    f"-> that is the n for the cost-of-matching result")
picked["in_both_arms"] = picked.dataset.isin(both)
```

The two arms do not contain the same datasets: matching is a search and succeeds to different
degrees, so a dataset can clear `min_pairs` in one arm and miss it in the other. Recording this
at selection time is why the R1 sample size is never a surprise later — this is the fix for the
95/187/189 confusion at its root.

---

# Part 5. `scripts/cloud_variants.py` — stage 11

## 5.1 The design, stated in the file

Stage in, run the existing code **unchanged**, stage out.

```python
def run_existing(what):
    cmd = [sys.executable, str(ROOT / "scripts" / "rehearsal_variants.py"), "--what", what]
    rc = subprocess.run(cmd, cwd=str(ROOT)).returncode
    if rc != 0:
        raise SystemExit(f"rehearsal_variants.py --what {what} failed (rc={rc})")
```

Reimplementing the assignment and conservation logic against a cloud filesystem would mean a
second copy of subtle code — strand handling, window offsets, ref-allele checks — diverging
silently from the copy that produced the published numbers. **Copying bytes is cheap; copying
logic is how two versions of a result appear.**

## 5.2 Staging only what is needed

```python
for r in panel.itertuples():
    dest = PROC_DIR / r.cell_line / r.protein / "dataset.tsv"
    if dest.exists():
        continue
    b = bucket.blob(f"processed/dinuc/{r.cell_line}/{r.protein}/dataset.tsv")
    if not b.exists():
        log(f"WARNING no processed dataset for {r.protein}:{r.cell_line}, skipping")
        continue
    b.download_to_filename(str(dest))
```

Only the study panel's datasets. Pulling all 244 would move data the stage never reads. The
`if dest.exists(): continue` makes a rerun cheap after a partial failure.

## 5.3 The marker discipline, again

```python
if len(sent) == len(OUTPUTS):
    bucket.blob(MARKER).upload_from_string(...)
else:
    log("NOT writing the completion marker: some outputs are missing, so a rerun "
        "must redo the work rather than skip it")
```

A marker written on partial success is worse than no marker: it makes the gap permanent and
silent.

---

# Part 6. `scripts/cloud_analysis.py` — stage 13

## 6.1 Differencing only the shared datasets

```python
m = gc.merge(dn, on="dataset", suffixes=("_gc", "_dn"))
m["cost"] = m.auroc_dn - m.auroc_gc
```

An **inner** merge, so only datasets present in both arms contribute. Differencing anything else
compares a dataset against nothing. This is `MATCHED = 187`, and it is why R1's n is 187 rather
than 189.

## 6.2 The ladder, computed once all three arms exist

```python
keep = set(arms["matched"].protein + ":" + arms["matched"].cell)
...
d = d.sort_values("delta", key=abs, ascending=False).drop_duplicates("vid")
d["block"] = (d.vid.str.split(":").str[0] + "_" +
              (d.vid.str.split(":").str[1].astype(int) // 1_000_000).astype(str))
f = cons.fit_delta_coef(..., blocks=d.block.to_numpy())
```

Line by line:

- `keep` restricts every arm to the datasets the matched arm covers, or the ladder compares
  **panels** rather than models.
- `drop_duplicates("vid")` after sorting by `|delta|` descending: a variant near several
  proteins' peaks appears several times, and those rows share a position, a label and a
  conservation value. Per-row inference treats ~19k independent observations as ~33k and reports
  intervals about a third too narrow. Keeping the largest-magnitude row per variant is the
  pre-registered collapse rule.
- `block` is a 1-Mb genomic bin, used as a gene proxy for the cluster bootstrap. Half the
  pathogenic variants sit in 7.7% of bins, so resampling rows rather than blocks understates
  uncertainty badly.

Conservation is added as its own rung with `coef = NaN`, because it has no delta coefficient by
construction. It is drawn as absent rather than zero.

---

# Part 7. `scripts/verify.py` — stage 14

## 7.1 Why tolerances, and what sets them

```yaml
r1_cost_of_matching:
  cost: {value: -0.1070, tol: 0.010}
```

Tolerances absorb what legitimately varies and nothing more:

- a different panel size (95 vs 187) shifts means by up to ~0.01
- BLAS thread count changes summation order and moves the last decimals
- GPU vs CPU inference differed by **max 1.14e-04 per variant**, measured on 19,051 shared
  variants with correlation 1.000000
- bootstrap intervals move with the seed in the third decimal

## 7.2 Checking claims, not just numbers

```python
at_least("gain RATIO (the paper's thesis)", gain_dn / gain_gc if gain_gc else None,
         spec["gain_ratio_min"])
```

The paper's thesis is that the gain over composition **at least doubles** under the harder
control. If that ratio inverts, the thesis is false regardless of how close every other number
lands. So the ratio is asserted directly, not left implicit in two means.

```python
ok = all(v is not None for v in vals) and all(a < b for a, b in zip(vals, vals[1:]))
record(ok, "model ordering holds", ...)
```

Ordering is the claim in R2 and R4, so it is checked as an ordering. `zip(vals, vals[1:])` is the
adjacent-pairs idiom.

```python
at_most("datasets significantly REVERSED", int((z < -1).sum()), spec["n_significantly_reversed_max"])
```

R3's strong form is "reversed on none". The **zero** is the part that matters, and it needs
per-dataset uncertainty — which is why `gini_sd` is carried through the cloud row. It was
originally dropped, and establishing that two apparent reversals were within 1 SE required a
separate local rerun purely because that number had been discarded on the way out.

## 7.3 The bug worth remembering

```python
# NOT `a or b`: pandas raises ValueError on DataFrame truthiness
d = T.get("matched_four_models.csv")
if d is None:
    d = T.get("matched95_four_models.csv")
```

`or` evaluates the left operand's truthiness and pandas refuses to define it for a DataFrame. The
verifier would have died with an unrelated error — which reads as "the verifier is broken" rather
than "the science did not reproduce".

---

# Part 8. `src/rbp/utils/cloud.py`

```python
def _resolve(explicit, env_key, cfg_key, what):
    if explicit: return explicit
    if os.environ.get(env_key): return os.environ[env_key]
    v = _from_config().get(cfg_key)
    if v: return v
    raise RuntimeError(f"{what} is not configured. ...")
```

Most-specific first: argument, environment, config. **No default**, because a default is how you
write into someone else's bucket, or read a stale one and never notice the run did nothing.

```python
@lru_cache(maxsize=1)
def _from_config(): ...
```

Called in loops; reading YAML per call is wasteful. `maxsize=1` because there is one config.

```python
def study_panel(bucket=None):
    b = (bucket or globals()["bucket"]()).blob(STUDY_PANEL_KEY)
    if not b.exists():
        return None
    ...
```

Returns **None**, not an exception, when no panel exists. Deliberate: prep itself runs before
selection and must process everything, and should not need to know it is special. `in_study_panel(None, ...)`
is therefore `True`.

---

# Part 9. `config/params.yaml` — the parameters that matter

| key | value | why this value |
|---|---|---|
| `cv.k` | 5 | five folds; `fold_roles` makes fold *f* test and *(f+1)%k* validation, so each fold is test once, val once, train k−2 times |
| `encode.min_pairs` | 400 | below this, a per-dataset AUROC interval is too wide to say anything. This threshold is why `MATCHED` is 187 and not 189 |
| `windows.size` | 101 | odd, so there is a true centre nucleotide for a variant to sit on |
| `variants.shifts` | [-40,-20,0,20,40] | the binding site's offset from the peak midpoint is unknown per variant, so five registers are scored |
| `variants.delta` | `max_abs` | maximum over shifts: "is there ANY register in which this variant matters". Averaging would dilute a real effect with windows that do not contain the site |
| `cloud.project_id` | `rbp-repro-2026` | the single value that redirects everything |

---

# Part 10. The test suite, and what each group is actually for

| file | pins |
|---|---|
| `test_no_hardcoded_project.py` | no source file names the author's project. Parametrised **per file**, so a failure names the offender. Searches `scripts`, `src`, `cloud`, `docker` — `docker` was missing, which is how two hardcoded ids survived a test written to prevent that |
| `test_panel_counts.py` | the four panel sizes and their nesting, plus that `DEEP` spans `FULL`'s size range so it cannot be a size-threshold sample |
| `test_locality_ism.py` | the probe on **constructed** signals where the answer is known: pure composition must give Gini ≈ 0, an implanted motif ≈ 0.92. The previous probe validated at r=+0.96 against a literature control and was still invalid, because both were driven by how strong the signal was |
| `test_train_folds.py` | fold role assignment. **Excluded from the CPU image**: it imports `rbp.train.data`, which imports torch |

The lesson from `test_locality_ism.py` generalises: **agreement with ground truth on a small
sample proves nothing if a third variable drives both.** Only constructed cases where the answer
is known by construction can validate an instrument.
