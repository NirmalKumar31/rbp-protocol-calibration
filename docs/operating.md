# 56. Operating and monitoring

How to actually run this: what to watch, how to read a failure, how to estimate cost before
spending it, and when to stop. Everything here is a technique that was needed at least once on
this project, with the incident that taught it.

---

# Chapter 1. The hierarchy of trust

When two sources disagree about what happened, believe them in this order. Every inversion of
this order on this project produced a wrong conclusion.

| rank | source | why | how it lies |
|---|---|---|---|
| 1 | **the artefact in GCS** | it either exists or it does not | it does not lie. It can be *incomplete* |
| 2 | **the completion marker** | written last, after the payload | absent marker + present payload = died between writes. This is the designed failure |
| 3 | **the task's own log lines** | what the code said it did | buffered output can lag, or vanish on SIGKILL |
| 4 | **Batch task counts** | the control plane's tally | lags reality by up to a minute; transitions look like regressions |
| 5 | **job state (RUNNING/FAILED)** | coarse | **FAILED can mean every real task succeeded** (bug 8) |
| 6 | **my own monitoring script** | derived | this lied three times on the original project |
| 7 | **the billing console** | authoritative on money | lags by hours. Useless during a run |

## 1.1 The three times monitoring lied, and what each teaches

**A zero-byte log.** A background job wrote nothing for twenty minutes and looked hung. Python
buffers stdout when it is not a terminal. `python -u` (or `PYTHONUNBUFFERED=1`, which the images
set) fixes it. **Teach: an empty log is not evidence of an idle process.**

**"Complete" at 273 of 475.** A wait loop grepped Modal's JSON for a completion field. The grep
failed to match, the loop treated no-match as done, and reported success at 57%. The fix is to
**assume alive on a parse failure** — `submit.sh` does exactly this with its `"")` case. **Teach:
a monitor's failure mode must be "keep waiting", never "declare victory".**

**78 tasks/hour when the real rate was 246.** A rate computed from two snapshots that were both
stale. **Teach: timestamp your samples and compute the delta over a measured interval**, which is
why every rate in this document comes from `T0`/`T1` arithmetic rather than a remembered earlier
number.

---

# Chapter 2. Watching a run

## 2.1 The one-liner for progress

Count **artefacts**, not task states. Artefacts are rank 1; task counts are rank 4.

```bash
gsutil ls "gs://${DERIVED}/processed/*/*/*/dataset.tsv" 2>/dev/null | wc -l
```

## 2.2 Measuring a rate honestly

```bash
A=$(gsutil ls "gs://${DERIVED}/processed/*/*/*/dataset.tsv" 2>/dev/null | wc -l | tr -d ' ')
T0=$(date +%s); sleep 300
B=$(gsutil ls "gs://${DERIVED}/processed/*/*/*/dataset.tsv" 2>/dev/null | wc -l | tr -d ' ')
T1=$(date +%s)
python3 -c "
import datetime
a,b,t=$A,$B,$T1-$T0
r=(b-a)/(t/60); left=488-b
print(f'{a} -> {b} in {t}s   rate {r:.2f}/min')
print(f'ETA {left/r:.0f} min -> ~{datetime.datetime.now()+datetime.timedelta(minutes=left/r):%H:%M}')
"
```

**Five minutes, not one.** Batch hands out tasks in waves; a one-minute window can show zero or
double the true rate.

**Never extrapolate from a single task.** I told the user prep would take 45–50 minutes based on
one log line showing a 47.5s task. It took 3.5 hours. The manifest is deliberately sorted
**biggest-first**, so early tasks are the *slowest* — throughput at the start is unrepresentative
by design. Measured over a window: 2.28/min.

## 2.3 Job state and counts

```bash
gcloud batch jobs describe "$JOB" --project="$PROJECT" --location="$REGION" \
  --format="value(status.state,status.taskGroups.group0.counts)"
# RUNNING   PENDING=416;RUNNING=6;SUCCEEDED=66
```

Counts that appear to go **backwards** are normal: a node cycling out returns its tasks to
PENDING. I once read that as a regression and invented a cause for it.

## 2.4 Events — read these first on any failure

```bash
gcloud batch jobs describe "$JOB" --format="json(status.statusEvents)" | python3 -c "
import json,sys
for e in json.load(sys.stdin)['status'].get('statusEvents',[]):
    print(e.get('type'), e.get('description','')[:150])
"
```

This is the rung most often skipped and the one that mattered most here.
`CODE_GCE_QUOTA_EXCEEDED` sat in this output while I theorised about spot preemption, switched
every job to on-demand, and measured no improvement.

`OPERATIONAL_INFO` is not an error — it is Batch telling you about the *infrastructure*, and it
is where quota, capacity and preemption all appear.

## 2.5 Task logs without drowning

```bash
gcloud logging read 'labels.job_uid="<job_uid>"' --project="$PROJECT" \
  --limit=20 --format="value(textPayload)" --freshness=25m | grep -vE "^report agent state"
```

**The filter is essential.** Batch agents emit enormous protobuf state dumps. An unfiltered read
of twelve entries returned several thousand tokens of agent metadata with two useful lines buried
in it. Filter to your own log format — every script here prints `[HH:MM:SS] message`, so
`grep -E "^\[[0-9]{2}:"` isolates them.

`--freshness` is required for recent logs; without it the default window can miss them.

## 2.6 What a healthy prep task looks like

```
[16:08:42] task 199: GEMIN5 K562 gc-matched  (61 KB of peaks)
[16:08:47] published gs://.../processed/gc/K562/SAFB2/  peak rss 0.17 GB
[16:08:47]   3744 pairs, 7488 rows in 47.5s
```

Three things worth noting. `peak rss 0.17 GB` against a 3500 MiB allocation means memory is
massively over-provisioned and could be cut to fit more tasks per node — except vCPU is the
binding constraint, so it would buy nothing. `3744 pairs` is checkable against the original run,
and matching it is the strongest single signal that the reproduction is real. And the task
announces *what* it is doing before doing it, so a hang tells you which dataset hung.

---

# Chapter 3. Reading a failure

## 3.1 The ladder, in cost order

Do not skip a rung. Each is cheaper than the next.

1. **Did the control plane accept it?** `gcloud batch jobs describe`. Rejected specs fail here
   with a precise field path — `Unknown name "networkInterfaces" at 'job.allocation_policy'`
   told me exactly what was wrong.
2. **What do the Events say?** §2.4.
3. **Are there VMs?** `gcloud compute instances list`. RUNNING job + no VMs = allocation
   failing.
4. **What is quota doing?** `gcloud compute project-info describe`. Compare **usage** to
   **limit**. Never filter to `limit > 0`: that is how a survey hid `GPUS_ALL_REGIONS = 0` from
   its own output, twice.
5. **Did the container start?** Look for `Runnable command line:` in the logs. Present means the
   image pulled and docker ran.
6. **What did the code say?** Your own log lines.
7. **Did it write anything?** `gsutil ls` the output prefix, and compare against the marker.
8. **Is the number right?** `scripts/verify.py`. Everything above can pass while the science is
   wrong.

## 3.2 Signature table

| symptom | almost always | first command |
|---|---|---|
| a `list` returns nothing, exit 0 | no project configured | `gcloud config get-value project` |
| job submits, dies in ~2 min | API not enabled / missing service agent | `gcloud services list --enabled` |
| RUNNING, no VMs, quota in Events | asking for exactly the quota limit | `project-info describe` |
| container crash-loops in seconds, 0 inputs | import error | read the container log |
| task hangs then hits `maxRunDuration` | outbound to a non-Google host, no route | check `EXTERNAL` for that stage |
| **FAILED job, complete output** | task count past the end of the manifest | compare `taskCount` to manifest rows |
| counts go backwards | a node cycled out; normal | nothing |
| budget shows $0 spent | `INCLUDE_ALL_CREDITS` | `gcloud billing budgets list` |
| 403 right after granting IAM | policy propagation | wait 30s, retry |
| `repository name must be lowercase` | an unexpanded `$VAR` reached docker | check substitutions |
| suite exits 0, too few tests | a collection error dropped whole files | run with `--collect-only` |
| `ValueError: truth value of a DataFrame` | `or` used on a DataFrame | explicit `is None` |

## 3.3 The two failures that look like success

These are the dangerous ones, because nothing turns red.

**A FAILED job whose work is complete.** `COUNT=189` against a 187-row manifest dispatched two
out-of-range tasks. Every real task succeeded; the job is red. **Check: compare `taskCount`
against the manifest length before believing a failure.**

**A SUCCEEDED job that did nothing.** 945 tasks reported "already complete, nothing to do"
because the image had not been rebuilt after a code change, so the completion-marker check saw
old markers. **Check: after any code change, confirm the digest changed.**

```bash
gsutil cat "gs://${PROJECT}-artifacts/images/cpu_digest.txt"
```

---

# Chapter 4. Cost

## 4.1 Estimate before spending, always

The probe pattern, from `cloud/modal/modal_variants.py`:

```python
@app.local_entrypoint()
def probe(index: int = 0):
    t0 = time.time(); rc = task.remote(index, True); el = time.time() - t0
    print(f"rc={rc} in {el:.0f}s")
    print(f"projected {N_TASKS} / {MAX_CONTAINERS}: "
          f"{el*N_TASKS/MAX_CONTAINERS/60:.1f} min wall, ${el*N_TASKS/3600*0.59:.2f} at T4 rates")
```

**One task, then multiply.** This caught all three Modal bugs — a missing module, a missing
dependency, and a 403 — for about a cent each, before a 94-task sweep.

Its estimate was also usefully **pessimistic**: it projected 8.9 minutes and $0.87; the real
sweep took ~2.5 minutes because containers stayed warm. Over-estimating from a cold start is the
right direction to be wrong in.

## 4.2 Eleven estimates, and the two that were badly wrong

Recorded because the *pattern* of error matters more than the numbers.

| guessed | measured | lesson |
|---|---|---|
| GPU speedup 100–200× | **29.6×** | never guess a hardware ratio |
| A100 = 6× T4 | **2.89×** | a 20M-parameter model does not saturate an A100 |
| bigger batch will help | **zero gain** | the bottleneck was elsewhere |
| RNABERT penalty = CNN's 1.65× | **4.9×** | compute-bound and memory-bandwidth-bound scale differently |
| Modal rate $1.10 + CPU + memory | **$1.10** | **44% too high.** GPU functions bundle CPU/RAM; they are not billed on top |
| prep in 45–50 min | **~3.5 h** | extrapolated from one task in a biggest-first manifest |

The Modal rate error is the instructive one. I built it by summing three published prices, which
*sounds* rigorous. The user caught it: "on the modal site it is saying i only spent 1.5$ extra
apart from 30$ and you are saying $36?" Back-solving from the dashboard gave the GPU price
alone. **When an estimate and a bill disagree, the bill is right.**

## 4.3 Fitting a cost model from finished work

Once some tasks have completed, stop estimating and measure:

```bash
# GPU-hours actually recorded, from the metrics each task uploads
gsutil cat "gs://${DERIVED}/runs/**/metrics.json" 2>/dev/null \
  | python3 -c "
import sys, json
tot=0
for line in sys.stdin:
    try: tot += json.loads(line).get('elapsed', 0)
    except Exception: pass
print(f'{tot/3600:.2f} recorded GPU-hours')"
```

`cloud/modal/guard.py` reports **two** numbers deliberately:

- an **upper bound**: elapsed wall time × max containers × rate
- a **lower bound**: recorded task time summed

The truth is between them. Reporting one number invites false precision — and at one point the
work-based estimate exceeded its own upper bound ($14.75 against $13.64), because per-task
overhead was counted for every task when Modal actually reuses warm containers.

## 4.4 The brakes, in order of how fast they act

| brake | latency | what it stops |
|---|---|---|
| `parallelism` / `max_containers` | **immediate** | the burn *rate*, before anything starts |
| `maxRunDuration` per task | minutes | a hung task billing forever |
| quota | immediate | absolute concurrency |
| a `confirm` gate | immediate | you, from starting it |
| budget alert | **hours** | nothing, on its own |
| billing killswitch | hours after the alert | everything, bluntly |

**A budget is a smoke alarm, not a sprinkler.** The brakes that matter act at submission time.

## 4.5 What this project actually cost

| | |
|---|---|
| GCP, original study | ~$8.55 of $300 credit |
| GCP, the rebuild so far | ~$3 of credit |
| Modal, SpliceBERT sweep | ~$31.50 (exhausted the $30 credit) |
| Modal, everything else | ~$1.40 out of pocket |
| **Real money, total** | **~$7.40** |

The whole study, including a full from-scratch reproduction, for under ten dollars of real money.
That is a design outcome, not luck: every expensive stage is preceded by a one-task probe and a
confirmation gate.

---

# Chapter 5. Deciding when to stop

## 5.1 Stop and investigate

- **Any destroy in a Terraform plan on a fresh project.** Additive by definition; a destroy means
  the state describes something else. This nearly cost the original study.
- **A rate that is an order of magnitude off** the probe's projection.
- **A SUCCEEDED stage that produced no artefacts.** Check the digest.
- **Any number outside `golden.yaml`'s tolerance.** The pipeline ran; the science did not
  reproduce. Do not write it up.

## 5.2 Do not stop

- **`OPERATIONAL_INFO: CODE_GCE_QUOTA_EXCEEDED`** when the job is progressing. Batch is retrying a
  node it cannot have. Annoying, harmless — and I spent real time on this.
- **Task counts going backwards.** A node cycled out.
- **A scatter of preempted tasks on spot.** Resumable; that is what markers are for.
- **A FAILED job with complete output.** Check the manifest length first.

## 5.3 The stop-loss discipline

Before starting a paid stage, write down what you expect. Afterwards, compare. If they differ by
more than 2×, stop and find out why before spending more.

This is not caution for its own sake. The 44% Modal rate error survived because nobody compared
the estimate to the bill until $6 of out-of-pocket spend had accumulated.

---

# Chapter 6. Runbooks

## 6.1 Resume after any interruption

```bash
cd rbp-repro
export GOOGLE_CLOUD_PROJECT=rbp-repro-2026
export PATH="../rna-binding-proteins/.venv/bin:$PATH"      # modal must be ON PATH
export PY=../rna-binding-proteins/.venv/bin/python

./run.sh status          # which artefacts exist
./run.sh from 6          # continue from the first unfinished stage
```

Every stage is idempotent. Rerunning a finished one costs seconds because markers are checked
first.

## 6.2 A stage failed; what now

```bash
# 1. which job?
gcloud batch jobs list --project=$GOOGLE_CLOUD_PROJECT --location=us-central1 \
  --format="table(name.basename(),status.state)" | head

# 2. events first
gcloud batch jobs describe <JOB> --format="value(status.statusEvents)" | tr ';' '\n' | tail -5

# 3. how much did it actually finish?
gsutil ls "gs://${DERIVED}/<prefix>/**" | wc -l

# 4. rerun. Completion markers skip what is done.
./run.sh stage N
```

## 6.3 Everything is on fire; stop the spend

```bash
# Batch
for j in $(gcloud batch jobs list --project=$P --location=us-central1 \
           --filter="status.state=RUNNING" --format="value(name.basename())"); do
  gcloud batch jobs delete "$j" --project=$P --location=us-central1 --quiet
done

# Modal
modal app list | grep rbp
modal app stop <APP_ID> --yes

# the blunt instrument, last resort
gcloud billing projects unlink $P
# reverse with:
gcloud billing projects link $P --billing-account=$BILLING_ACCOUNT
```

Then verify nothing survived:

```bash
gcloud compute instances list --project=$P            # should be empty
modal app list | grep -v stopped                       # should be empty
pgrep -f "modal run|run.sh"                            # should be empty
```

## 6.4 Verify a reproduction

```bash
./run.sh stage 14                       # reads tables from GCS
python scripts/verify.py --local results/tables    # or from a local directory
```

Exit 0 means every golden number reproduced. Exit 1 names the claim that broke.

## 6.5 Rotate the Modal credentials

Do this when a run finishes, and immediately if a token was ever pasted into a chat or a log.

```bash
# modal.com/settings/tokens -> revoke, create new
modal token set --token-id ak-NEW --token-secret as-NEW --profile=NAME

# and the GCP key inside the secret
gcloud iam service-accounts keys list \
  --iam-account=rbp-modal@$GOOGLE_CLOUD_PROJECT.iam.gserviceaccount.com
gcloud iam service-accounts keys delete <OLD_KEY_ID> \
  --iam-account=rbp-modal@$GOOGLE_CLOUD_PROJECT.iam.gserviceaccount.com
```

The `rbp-modal` key is the **only credential that leaves Google's network**, which is why its IAM
condition restricts it to three prefixes. Rotating it is cheap; assuming it is safe is not.

---

# Chapter 7. The operating rules, condensed

1. **Trust artefacts over status.** An object either exists or it does not.
2. **Measure, never extrapolate.** One sample is not a rate.
3. **Read the Events before theorising.** The cause is usually already written down.
4. **A monitor must fail toward "keep waiting"**, never toward "done".
5. **Probe one task before fanning out.** A cent to learn what a dollar would teach.
6. **Compare the bill to the estimate.** The bill is right.
7. **Cap concurrency, not just budget.** Budgets lag by hours; parallelism acts now.
8. **Never type a count.** Read it from the manifest.
9. **Write the marker last.** Make the survivable failure the one that happens.
10. **Verification is a stage.** A pipeline that completes and produces different science is
    worse than one that crashes.
