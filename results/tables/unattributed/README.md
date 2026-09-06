# Tables with no producing script

Two files, moved here rather than deleted or left in place.

`scripts/provenance.py` walks every committed table and finds the script that writes it. Three
had none. One, `positive_set_overlap.csv`, is quoted in the Discussion, so its producer was
reconstructed and its measure turned out to be wrong; see `scripts/positive_set_overlap.py`.

These two are the others. Neither is cited in the manuscript, in `SUBMISSION.md`, in
`config/golden.yaml` or in `scripts/verify.py`, and no script in the repository produces them.

| file | what it appears to be |
|---|---|
| `backend_replication.csv` | 20 datasets, GC and dinucleotide contributions recomputed against the published value. A BLAS-backend or platform replication spot check |
| `substitution_baseline.csv` | 84 datasets, `auroc_subtype` and `auroc_kmer_delta`. ClinVar variant scoring, which belongs to the earlier study rather than to this paper |

They are kept because they are evidence of work that was done and deleting evidence is not an
improvement. They are moved out of `results/tables/` because everything in that directory is
supposed to be regenerable by a named command, and a reader cannot tell by looking which files
carry that guarantee. Anything here does not.
