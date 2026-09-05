# What this cost, and what a rerun would cost

One place, because there were five and they disagreed. An external review found the release
quoting roughly $24, $34, $37 to $38, $50 and under $11 for what a reader would reasonably read
as the same quantity. They were not the same quantity: some were the original variant-scoring
study, some this paper's three-arm sweeps, some a budget ceiling, some a burn rate, and one was
a forecast made before the GPU rate was corrected by 44%. None of them said which.

Four different things, kept apart:

## 1. Measured, the published run

What was actually billed for the work behind this paper. Sources named per row; nothing here is
a forecast.

| component | measured | source |
|---|---:|---|
| GCP, preprocessing and ingest, original study | ~$8.55 | billing console, against $300 trial credit |
| GCP, the three-arm rebuild | ~$3 | billing console, against the same credit |
| Modal, the first SpliceBERT sweep | ~$31.50 | Modal dashboard; exhausted a $30 credit |
| Modal, the three-arm × three-model sweeps | $18.97 | Modal dashboard, against a $18.91 forecast |
| **Real money out of pocket** | **~$20** | the rest was trial and platform credit |

The two GPU rates the forecast used, back-solved from recorded A10G seconds in `metrics.json`
across the bias-aware arm's 940 runs: **$14.58 per million pairs for the CNN** and **$27.25 per
million pairs for SpliceBERT**, both at the A10G list price of **$1.10/GPU-hour**. Those are the
numbers `cloud/modal/modal_gc_sweep.py` prices a manifest with, and the $18.97 against $18.91
above is how well they held.

`results/tables/device_portability.csv` carries one of these as a committed assertion: the CNN
sweep at $6.73 over 470 fold-runs.

## 2. Forecast, a rerun today

Not the same as row 1, because row 1 includes a sweep that was later superseded and excludes
work that was done under credit before the rate was corrected.

| stage | where | forecast |
|---|---|---:|
| images, ingest, panel, preprocess | GCP Batch | ~$3 |
| CNN, three arms | Modal A10G | ~$20 |
| SpliceBERT, three arms | Modal A10G | ~$37 |
| **total** | | **~$60** |

Higher than what was spent, and that is the honest direction: the published run got a $30 Modal
credit and a $300 GCP trial credit that a reader reproducing it will not have.

## 3. Budgets and kill thresholds

These are ceilings, not estimates, and they are deliberately well above the forecast so that
ordinary estimation drift cannot trip them. A guard that fires on a 20% miss destroys a
partially complete sweep for no benefit; see the note in `cloud/modal/guard.py`.

| control | value | what it does |
|---|---:|---|
| `cloud/modal/guard.py --budget` | $40 | stops the Modal app on the upper-bound estimate |
| Modal burn at full fan-out | $11.00/h | 10 A10G containers at $1.10 |
| GCP budget alert | $30 | email only, stops nothing |
| `cloud/killswitch` threshold | $40 | detaches the project's billing account |

## 4. Nothing above is required to check the paper

`python scripts/verify.py --local results/tables` runs the 990 published assertions against
committed tables. It costs nothing, needs no account, and is the only thing a reader has to run.
