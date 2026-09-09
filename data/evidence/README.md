# What is in here

Per-example out-of-fold model scores: the evidence that makes the model-class comparison
recomputable from this repository alone, rather than asserted against a summary table.
`scripts/recompute.py` rebuilds 285 published AUROCs from these files, and
`results/tables/PROVENANCE.csv` records which released table each one supports.

No genomic sequence is redistributed. A row is a window identifier, a label, a fold and a
score; see `LICENSE` beside this file.

## The three protocol arms

One directory per negative-set protocol, identical layout:

    <arm>/<cell>/<protein>/<model>/fold<0-4>/scores.tsv.gz

| directory | protocol | datasets | files |
|---|---|---:|---:|
| `scores/` | dinucleotide-matched | 95 | 950 |
| `scores_gc/` | GC-matched | 94 | 940 |
| `scores_neg2/` | bias-aware, other RBPs' sites | 94 | 940 |

`models` are `cnn` and `splicebert`, five folds each, so a dataset contributes ten files per
arm. **`scores/` carries 95 datasets and the other two carry 94.** That is the panel difference
the paper reports throughout, not a gap: NCBP2 in K562 clears the minimum-pair floor under
dinucleotide matching and not under GC matching, so it has no three-protocol result.
`docs/PANELS.md` is the single place that reconciles every count.

Columns: `id`, `label`, `fold`, `score`. `id` is the window identifier, `label` is 1 for a
positive and 0 for its matched negative, `fold` is the frozen chromosome-blocked assignment,
and `score` is the out-of-fold prediction. The 4-mer is not here because it is refit from the
window tables rather than scored once.

## The other two directories

`rehearsal/<cell>/<protein>.scores.tsv.gz`, 95 files, same four columns. The k-mer rehearsal
pass over the whole panel that ran before the neural sweeps. `scores/rehearsal/<cell>/<protein>.json`,
95 files, is that pass's per-dataset summary, one JSON object with the AUROC and its interval.

`scores_md/`, 474 files named `<cell>_<target>__<cell>_<donor>.csv`. The multi-donor
wrong-protein control described in `scripts/multidonor_analysis.py`, five donors per target.
Columns: `protein, cell, vid, label, fold, delta, platform, accelerator, weights_from`. This
belongs to the earlier variant-scoring study and no claim in the paper depends on it; README.md
at the repository root explains why that study is kept rather than deleted.
