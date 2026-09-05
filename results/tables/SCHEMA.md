# What the columns in `results/tables/` mean

124 CSVs, two shapes. An external review asked for a data dictionary on the grounds that an
archival reader should not have to read the producing script to learn what a column holds, which
is right, and the answer is short because the shapes are conventional rather than per-file.

## Shape 1: a summary table, one row per asserted quantity

44 of the tables. `scripts/verify.py` reads these, and `scripts/audit_manuscript.py` treats every
cell of a table with a `check` column as a value the manuscript is allowed to quote.

| column | type | meaning | missing |
|---|---|---|---|
| `check` | string | the quantity, in words. The primary key: unique within a file, and the string `verify.py` looks the row up by, so renaming one is a breaking change | never |
| `value` | float | the quantity | never |
| `ci_low`, `ci_high` | float | 95% interval bounds. **Protein-clustered percentile bootstrap, 4000 draws**, unless the `note` says otherwise: the 15 proteins measured in both cell lines contribute two correlated rows and are resampled together | empty where the quantity is a count, a ratio of means, or otherwise has no resampling distribution |
| `n` | int | datasets, or rows, behind the value. Usually 94 | never |
| `note` | string | scope, caveat, or the definition of a non-obvious quantity | empty where there is nothing to say |

Column order varies between files and carries no meaning; read by name.

## Shape 2: a per-dataset table, one row per dataset

Named `*_per_dataset.csv`, always beside the summary table that aggregates it. These are the
evidence: the summary is recomputable from them, which is what `--from-cache` does.

| column | type | meaning |
|---|---|---|
| `dataset` | string | `PROTEIN:CELL`. The primary key |
| `protein` | string | ENCODE target. **Not unique** — 15 appear in both cell lines, which is why intervals are protein-clustered |
| `cell` | string | `K562` or `HepG2` |
| `n_{arm}` | int | rows in that arm's dataset, which is twice the pair count |
| `comp_{arm}` | float | composition-only pooled out-of-fold AUROC in that arm |
| `gain_{arm}`, `delta_{arm}` | float | nested contribution: AUROC(19 composition features + model score) − AUROC(composition alone), both pooled out of fold |

## The arm suffixes

One protocol per suffix, used in column names throughout.

| suffix | protocol |
|---|---|
| `gc` | GC-matched negatives |
| `dn`, `dinuc` | dinucleotide-matched negatives |
| `neg2` | bias-aware: negatives are other RBPs' binding sites |
| `neg2_rm` | bias-aware, donor draw stratified on transcript region |
| `shuf` | dinucleotide-shuffled |

Both `dn` and `dinuc` occur. `dn` is the column-name form and `dinuc` the directory and
`--arm` form; they are the same protocol. The review flagged the drift and it is real, but
renaming a column breaks every `check` string in `golden.yaml`, so it is documented here rather
than churned.

## Units and conventions

- Every AUROC is **pooled** over the five out-of-fold score vectors, not averaged per fold. See
  `scripts/auroc_aggregation.py`, which measures what that choice costs.
- Contributions are AUROC differences, so they are dimensionless and on `[-1, 1]`.
- Costs are US dollars; see `docs/COST.md`.
- An empty cell means not applicable. `NaN` means computed and undefined — a ratio whose
  denominator is not positive, or an AUROC on a fold with one class.
- Nothing here is rounded for presentation. The manuscript rounds; these do not.

## Provenance

Every table names its producing script in the `run.sh` stage that writes it, and
`scripts/verify.py` asserts its published values against `config/golden.yaml`. Two tables are
rebuilt end to end rather than compared against a record: `recompute.py` rebuilds 285 AUROCs
from committed per-window scores, and `k_sweep.py` rebuilds the headline contrast from raw
sequence.
