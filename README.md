# RBP composition confound: reproducible rebuild

Negative-set construction, not architecture, determines what an RBP binding benchmark
measures. Correcting it lowers headline scores while **doubling** the share attributable to
learning, and that learned signal is positionally localised and partly transfers to disease
variants.

| | result |
|---|---|
| **R1** | Dinucleotide-matched negatives cost 0.107 AUROC vs GC-matched, on every dataset. Gain over a composition-only baseline **rises** +0.027 to +0.064 |
| **R2** | composition 0.628 < k-mer 0.688 < CNN 0.708 < SpliceBERT 0.809, identical splits |
| **R3** | SpliceBERT's sensitivity is more positionally concentrated; reversed on zero datasets |
| **R4** | ClinVar ladder: k-mer 0.529 < wrong-protein head 0.656 < right-protein head 0.829 (conservation 0.906) |

## Run it

```bash
export GOOGLE_CLOUD_PROJECT=your-new-project
./run.sh preflight     # spends nothing, gates everything
./run.sh all           # pauses at each paid stage
./run.sh stage 14      # verify against config/golden.yaml
```

Full procedure: **[docs/REPRODUCE.md](docs/REPRODUCE.md)**.
Why the dataset counts differ: **[docs/PANELS.md](docs/PANELS.md)**.

## Cost

~$5 GCP credit, ~$32 Modal. Stage 9 (SpliceBERT) is 95% of it. Every paid stage asks first.

## Design rules

1. **No local compute.** Every stage runs in a container on Batch or on Modal.
2. **No hardcoded project id.** Everything resolves through `rbp.utils.cloud`; a test fails
   the build if a literal reappears.
3. **The panel is an artefact, not a flag** (`manifest/study_panel.tsv`), written once.
4. **Task counts come from manifests**, never typed.
5. **Completion markers are written last**, so an interrupted stage redoes its work.
6. **Verification is a stage.** Reproducibility that is not checked is not reproducibility.
