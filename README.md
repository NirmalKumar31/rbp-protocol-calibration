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
| **R1** | Three protocols, one model: nested contribution **+0.0663** (dinucleotide-matched), **+0.0265** (GC-matched), **+0.0122** (negatives = other RBPs' sites). The model's own AUROC moves the OPPOSITE way, falling from 0.7981 to 0.6879 in 94/94 while the contribution rises 2.5x |
| **R1m** | No monotone rescaling reaches protocol independence; floor **2.00x**, achieved by dividing by the baseline's own headroom |
| **R1l / R1n** | Protocol and baseline are largely confounded, but not unmeasurably so. Given the baseline, the protocol label adds **1.0%** of variance; given the protocol, the baseline adds **11.0%**. The 0.0056-wide common support was the THREE-WAY intersection; pairwise, GC vs neg2 overlaps **0.212** over **130/188** cells, and there a protocol residual survives at matched baseline, **−0.0081** [−0.0130, −0.0036] |
| **R1g** | Holds for three model classes: k-mer **+0.0398**, CNN **+0.0530**, SpliceBERT **+0.0864**. The multiplier is mostly a property of the PROTEIN, and the direct test is the evidence: the same protein's log multiplier agrees across cell lines at **r = +0.586** (p = 0.0001, 40 protein x model pairs). As a variance share it is 64.8%, but protein is nearly the *dataset* factor (68.0%), so against a null that permutes protein labels between datasets the excess is **+7.2** points (p = 0.0030), not the +35.1 a wholesale permutation suggests. Cell line is not detectable (15 informative proteins); model class is small but **not** null (2.8%, p = 0.023) |
| **R1o / R1r** | Raising the composition baseline to order 3 absorbs a near-CONSTANT absolute amount from every model class (**+0.021** GC, **+0.054** dinuc), which is nearly all of a 4-mer's contribution and a quarter of SpliceBERT's. Over a trinucleotide baseline the 4-mer is positive in only **65/94** datasets against SpliceBERT's **94/94**. The protocol contrast survives at order 3 for all three models |
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
