# The four panels, and why they are different sizes

Every count in this project is one of four panels. They are nested, each difference has one
cause, and no number is a mistake. This file is the single place that says so; anything that
disagrees with it is wrong.

| n | name | definition | why not larger |
|---|---|---|---|
| **189** | `FULL` | every dataset in the dinucleotide arm clearing `min_pairs: 400` | this is the whole panel |
| **187** | `MATCHED` | datasets present in **both** the gc and dinuc arms | `DDX51:K562` and `NCBP2:K562` matched only 360 and 384 pairs under GC matching, below the 400 floor |
| **95** | `DEEP` | datasets with all four models | `--every 2` systematic sample of `FULL`. A cost decision: SpliceBERT on 189 was unaffordable |
| **94** | `VARIANT` | `DEEP` datasets carrying ClinVar variants | `NCBP2:K562` has zero ClinVar variants near its peaks |

    VARIANT (94) ⊂ DEEP (95) ⊂ FULL (189)     MATCHED (187) ⊂ FULL (189)

## Which panel answers which question

| claim | panel | why that one |
|---|---|---|
| cost of proper matching, −0.1070 | `MATCHED` 187 | needs both arms; a dataset in only one cannot be differenced |
| four-model comparison | `DEEP` 95 | the only set where all four models exist on identical splits |
| k-mer / CNN / composition at scale | `FULL` 189 | those three ran everywhere |
| ISM locality | `DEEP` 95 | needs SpliceBERT weights |
| ClinVar | `VARIANT` 94 | needs both SpliceBERT weights and variants |

## The two rules that prevent the confusion

**1. Never compare across panels without saying so.** The gc-vs-dinuc difference must use the
187 both arms share, not 189 minus 187. The k-mer arm in the variant comparison is recomputed
on the 94, not quoted from the 187-dataset rehearsal, or the model would be confounded with
the panel.

**2. `DEEP` is a sample, so quantify the sampling.** Systematic by pair rank (`--every 2`),
which is unbiased in size by construction — a size-threshold sample would confound the subset
with dataset size, and AUROC correlates with size at r = +0.53 to +0.67. Measured: `DEEP` sits
within ~0.008 AUROC of `FULL` on every quantity, and the sign is NOT consistently optimistic
(k-mer −0.0010, composition −0.0023, CNN −0.0065, median CNN gain +0.0083). An earlier draft
quoted "+0.0082 optimistic" as a blanket caveat; that figure is one statistic, not the rule.

## Protein counts, which are also not a mistake

Datasets are protein x cell line, so protein counts are always lower: `FULL` 132 proteins,
`MATCHED` 131, `DEEP` 79. A protein assayed in both K562 and HepG2 is two datasets and one
protein, and the CV grouping is by chromosome, so the two are not independent replicates of
each other in any statistical sense -- they share the genome, not the assay.
