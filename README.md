# RBP negative-set calibration: a reproducible rebuild

**What a sequence model appears to contribute is set by how the benchmark's negatives were
built.** Holding the model, the positives, the folds and the estimator fixed and changing only
the negative-set protocol, the same 4-mer's contribution over a composition baseline moves
**5.4-fold** (95% CI 4.4 to 6.6). No rescaling of AUROC removes it: across eight monotone
transforms the range never falls below **2.00x** [1.67, 2.46].

**Report the composition-only AUROC under the same negative-set protocol alongside every
headline AUROC, and never compare contributions measured under different protocols.**

| | result |
|---|---|
| **R1** | Three protocols, one model: nested contribution **+0.0663** (dinucleotide-matched), **+0.0265** (GC-matched), **+0.0122** (negatives = other RBPs' sites). Apparent AUROC moves the OPPOSITE way, falling 0.1095 in 94/94 while the contribution rises 2.5x |
| **R1m** | No monotone rescaling reaches protocol independence; floor **2.00x**, achieved by dividing by the baseline's own headroom |
| **R1l** | Protocol and baseline are confounded by construction: their baseline ranges overlap in a window 0.0056 AUROC wide containing 3 of 282 cells |
| **R1g** | Holds for three model classes: k-mer **+0.0398**, CNN **+0.0530**, SpliceBERT **+0.0864**. The multiplier is mostly a property of the PROTEIN: **64.8%** of the log multiplier's variance, against the **29.7%** any 79-level factor absorbs from relabelled data (excess **+35.1**, p < 0.0005). Cell line is noise (0.2%, p = 0.49); model class is small but **not** null (2.8%, p = 0.023) |
| **R1r** | The order-3 collapse is the 4-mer's alone: surviving a trinucleotide baseline, k-mer **21.7%**, CNN **60.5%**, SpliceBERT **75.0%**, intervals non-overlapping |
| **R1c / R1j** | Two artifacts bounded against matched placebos: strand **−0.0055** (85% survives), untranscribed negatives **−0.0043** (90% survives) |

## Run it

```bash
export GOOGLE_CLOUD_PROJECT=your-new-project
./run.sh preflight     # spends nothing, gates everything
./run.sh all           # pauses at each paid stage
./run.sh stage 15      # verify against config/golden.yaml
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
