# The panels, and why they are different sizes

Every count in this project is one of a small number of panels. They are nested, each
difference has one cause, and no number is a mistake. This file is the single place that says
so; anything that disagrees with it is wrong.

**This file describes THIS pipeline.** An earlier, discarded build had a different structure
(`FULL` 189 datasets, `MATCHED` 187, losing `DDX51:K562` and `NCBP2:K562` to GC matching).
Those numbers are gone and must not reappear. The digit 189 survives into this build meaning
something completely different, which is exactly the kind of collision this file exists to
prevent — see the warning below.

## The panels

| n | name | definition | why not larger |
|---|---|---|---|
| **189** | `CANDIDATE` | every dataset in the dinucleotide arm clearing `min_pairs: 400` | 101 K562 + 88 HepG2; the pool the study panel was drawn from |
| **95** | `STUDY` | the analysed panel | systematic sample of `CANDIDATE` by pair rank. A cost decision: SpliceBERT on 189 was unaffordable |
| **94** | `GC` | `STUDY` datasets that also clear the floor under GC matching | `NCBP2:K562` matches 384 pairs under GC and 406 under dinucleotide, so it clears 400 in one arm only |
| **95** | `VARIANT` | `STUDY` datasets with ClinVar variants near their peaks | all of them, in this build |
| **82** | `VARIANT-USABLE` | variant datasets with both classes present, so a per-dataset AUROC exists | 13 carry only one label class; an AUROC does not exist there, and a zero is not a substitute |
| **44** | `VARIANT-POWERED` | ≥20 pathogenic variants | median coverage is 140 variants per dataset and many carry a handful of pathogenic ones, so their per-dataset AUROC is noise |

    GC (94) ⊂ STUDY (95) = VARIANT (95) ⊃ VARIANT-USABLE (82) ⊃ VARIANT-POWERED (44)
    STUDY (95) ⊂ CANDIDATE (189)

## THE 189 COLLISION, which has already caused one wrong assertion

**189 means two different things in this project and neither is a panel of analysed
datasets.**

1. `CANDIDATE`, the 189-dataset pool the study panel was sampled from.
2. The **stage 7 task count**: 95 dinucleotide tasks + 94 GC tasks, one per dataset per arm.

These are numerically equal by coincidence. A test in `tests/unit/test_panel_counts.py` once
asserted a 189-dataset analysed panel, inherited from the discarded build, and it survived
because the tests that would have contradicted it were skipping on a missing file. Both the
task count and the pool size are now asserted separately and explicitly.

## Which panel answers which question

| claim | panel | why that one |
|---|---|---|
| composition share, 94.8% → 67.8% | `GC` 94 | paired across arms; a dataset in only one cannot be differenced |
| cost of matching, −0.1095 | `GC` 94 | same reason |
| four-model comparison | `STUDY` 95 | the only set where all four models exist on identical splits |
| ISM locality | `STUDY` 95 | needs SpliceBERT weights |
| ClinVar specificity, +0.0645 | `VARIANT-POWERED` 44 | a per-dataset AUROC needs both classes and enough of the rare one |
| the trivial positional baseline | `VARIANT-POWERED` 44 | compared against the model on identical datasets |

## The rules that prevent the confusion

**1. Never compare across panels without saying so.** The GC-vs-dinucleotide difference uses
the 94 both arms share, not 95 minus 94.

**2. `STUDY` is a sample, so quantify the sampling.** Systematic by pair rank, taking every
second row of the size-sorted candidate list. That spans the full size range by construction —
a size-threshold sample would confound the subset with dataset size, and AUROC correlates with
size at r = +0.13 (composition) to +0.67 (SpliceBERT). Measured: `STUDY` spans `CANDIDATE` size
percentile **0–99**, shown in `f0`.

It is a deterministic cost-driven subset, not a probability sample, and the earlier wording
here ("unbiased in size by construction") overstated that. The phase is fixed at row 0 rather
than drawn, so the procedure guarantees range coverage but cannot be assumed free of aliasing
against any structure that happens to be ordered by pair count or by the tie-break. Every
interval in the paper is conditional on this panel.

**3. Report the unpowered stratum too.** The ClinVar specificity result is claimed on the 44
adequately powered datasets. The all-82 stratum shows nothing (gap −0.011, p=0.87) and is
asserted in `golden.yaml` precisely so that it is never quietly omitted. Stratifying by power
is legitimate here because the effect grows with power (rho=+0.52) while the wrong-protein
floor stays flat — but the reader gets to check that, which means printing both.

## Protein counts, which are also not a mistake

Datasets are protein × cell line, so protein counts are always lower: `STUDY` is 95 datasets
and **79 proteins**, of which **16** appear in both K562 and HepG2. A protein assayed in both
lines is two datasets and one protein, and the CV grouping is by chromosome, so the two are
not independent replicates of each other in any statistical sense — they share the genome,
not the assay.
