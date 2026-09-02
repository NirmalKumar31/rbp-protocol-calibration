# Materials and Methods

## Panel selection

Candidate experiments were every released ENCODE eCLIP experiment in K562 and HepG2 (139 and
105 respectively). Where a protein had several peak files we took the one with the most
biological replicates, breaking ties by ENCODE's `preferred_default` flag, and required
`output_type=peaks`, `assembly=GRCh38`, `status=released` and `file_format=bed`. Every candidate
was preprocessed under both composition-matching protocols; datasets retaining at least 400
out-of-fold scored pairs were eligible.

The study panel was then defined once by the dinucleotide arm, which is the stricter control:
eligible datasets were sorted by pair count with a stable sort and every second one retained,
giving 95 datasets over 79 proteins. Sampling systematically by rank rather than taking the
largest N matters because AUROC correlates with dataset size at r = 0.53 to 0.67 in these data,
so a size-truncated panel would confound the panel with the measured quantity; the retained
panel spans the 0th to 99th percentile of the candidate size distribution. Of the 95, 94 also
cleared the floor in the GC arm and constitute the paired panel used for every three-protocol
result, giving 282 protocol-dataset cells; the ladder in Table 1 uses all 95. The panel was
written once to a manifest and every downstream stage reads it, so membership is not decided
anywhere else. Supplementary Table S1 lists every dataset with its ENCFF file and ENCSR
experiment accession.

Fifteen proteins are assayed in both cell lines and contribute two datasets each. All headline
intervals therefore resample proteins rather than datasets (below).

## Windows

Positive windows are 101 nt centred on the midpoint of each peak interval. Peak width and all
significance columns are deliberately unused: no p-value, q-value or fold-change filter is
applied, so the positive set is every called peak. Minus-strand windows are reverse-complemented
and transcribed to RNA, so every sequence is the strand the protein sees. Windows containing an
ambiguous base, falling outside the assembly, duplicating an existing window by start
coordinate, or carrying no region annotation are dropped and counted. Analysis is restricted to
chr1 to chr22 and chrX; chrY and all scaffolds and alternates are excluded.

Region class is assigned from GENCODE v45 as the first of (5'UTR, 3'UTR, CDS, non-coding exon,
intron) whose interval set contains the window midpoint. Region intervals are built from all
transcripts with no `gene_type` filter, merged, and strand is deliberately discarded at this
step; introns are computed as the gaps between consecutive exons. GTF 1-based inclusive
coordinates are converted to 0-based half-open.

## The three negative-set protocols

Each protocol produces one negative per positive. All three use a single random stream seeded at
7 per protein.

**GC-matched.** Candidates are windows of the same region class on the same chromosome, at least
500 nt from any peak of the target protein, drawn with probability proportional to the length of
each free interval. A candidate is accepted if its GC content is within 0.05 of the positive's.
The matcher relaxes: after 25 unsuccessful draws the tolerance doubles to 0.10, and if 40 draws
fail the best candidate seen is accepted provided it is within three times the nominal
tolerance, that is 0.15. **The operative acceptance bound is therefore 0.15, not 0.05**, and we
report achieved quality below rather than the nominal figure. Only the target protein's own peaks
are excluded, so a GC-matched or dinucleotide-matched negative may be another protein's binding
site. The negative inherits the positive's strand label; the consequences are quantified in the
Results.

**Dinucleotide-matched.** Positives are grouped by region and chromosome, and for each group a
candidate pool of `max(1500, 8 x n)` windows is sampled. For each positive the nearest unused
candidate is assigned greedily under the L1 distance over the sixteen dinucleotide **counts**,
searched with a k-d tree over the forty nearest neighbours. Counts rather than frequencies are
used because integer distances are exactly representable, and each candidate additionally
carries a tiebreak column that increases strictly with genomic position and is bounded below
0.5, so it can never reorder candidates whose integer distances differ but resolves exact ties
by position. Without both devices the choice among tied candidates depends on the floating-point
behaviour of the k-d tree implementation and a different negative is selected for about 8% of
pairs between CPU architectures. Assignment is greedy rather than optimal; achieved distances
are reported so the quality of the approximation is visible. Used coordinates are tracked across
all groups, because one interval can be annotated as both a non-coding exon and a 3'UTR and so
appear in two region pools. No maximum distance is imposed.

**Bias-aware (other RBPs' sites).** Following [3], negatives for a target are the positive
windows of the other panel proteins assayed in the same cell line, sampled 1:1 without
replacement, excluding any donor window within 500 nt of any of the target's own positive
windows. Two deviations from [3] should be noted: our donor pool is the study panel (40 to 48
donor proteins per dataset, median 45) rather than a full ENCODE survey, and exclusion is by
distance from the target's 101 nt windows rather than from its full peak intervals. Sampling is
performed within fold, so every pair stays inside one cross-validation fold. These negatives are
transcribed, accessible to crosslinking and strand-correct by construction, and no composition
matcher touches them; they are not composition-matched in any respect, which is central to
interpreting the results.

## Achieved match quality

Because the paper argues from what the matchers control, we report what they achieved over all
456,734 pairs per arm rather than their specifications.

| arm | median abs. GC gap | p90 | p99 | max | within 0.05 | median dinucleotide L1 |
|---|---|---|---|---|---|---|
| GC-matched | 0.0297 | 0.0495 | 0.1089 | 0.1486 | 94.8% | 0.500 |
| dinucleotide-matched | 0.0198 | 0.0693 | 0.1485 | 0.4159 | 85.0% | 0.220 |
| bias-aware | 0.1387 | 0.3169 | 0.4555 | 0.6733 | 20.4% | 0.700 |

Dinucleotide L1 is in frequency units on a 0 to 2 scale. Two points follow. The GC matcher holds
its nominal tolerance for 94.8% of pairs and its observed maximum, 0.1486, sits at the 0.15
fallback cap. And "dinucleotide-matched" means a median L1 of 0.220 rather than zero: the
matcher improves composition distance 2.27-fold relative to the GC arm but does not eliminate
it, so statements about degrees of freedom controlled are exact for the design and approximate
for the realisation. The bias-aware arm is not composition-matched at all, with a median GC gap
4.7 times the GC arm's.

## Cross-validation

Evaluation is grouped 5-fold cross-validation over whole chromosomes. The chromosome-to-fold map
is frozen once for all datasets in both cell lines, so the fifteen proteins measured in both
lines are evaluated on the same chromosomes and a between-line difference cannot be a partition
artefact. The map was chosen by random-restart hill climbing on peak counts only, never on
labels or model output, minimising the summed squared deviation from equal mass per fold with
every dataset weighted equally; it was required to beat round-robin assignment, to hold the
median per-dataset maximum deviation at or below 0.05, and to place no more than half of any
dataset's data in one fold. Fold *i* is the test fold, fold *i*+1 modulo 5 is validation and the
remaining three train, so each fold is test once, validation once and training three times, and
no additional data is spent on model selection.

Every window is scored exactly once out of fold and all AUROCs are computed on the pooled
out-of-fold score vector. Pooling rather than averaging per-fold AUROCs makes the estimate
invariant to fold size, which matters because the least balanced dataset places 0.42 of its data
in a single fold.

Leakage between folds was audited by exact 32-mer sharing rather than by clustered identity: a
101 nt window contains 70 distinct 32-mers and two unrelated sequences share one with negligible
probability, so a shared 32-mer is evidence of common origin.

## The composition baseline and the nested contribution

The composition baseline is a 19-column design matrix: three mononucleotide frequencies, fifteen
dinucleotide frequencies and the Shannon entropy of the window's base composition. One level is
dropped from each frequency family because frequencies sum to one and retaining all of them makes
the design singular; dropping a level does not change the space the block spans. Every column is
standardised. GC content is not a feature: it is the spanned combination C+G, and because
mononucleotide counts are marginals of dinucleotide counts the block occupies the fifteen
degrees of freedom of the dinucleotide simplex. The GC-matched protocol therefore constrains one
linear functional of that space and the dinucleotide-matched protocol constrains all fifteen.
**The baseline is refit on each protocol's own windows and is never carried between protocols.**

The nested contribution of a model is

> AUROC(composition + standardised model score) minus AUROC(composition alone),

both fitted by L2-penalised logistic regression (C = 1, lbfgs, at most 3000 iterations) once per
fold on that fold's training rows, both evaluated as pooled out-of-fold linear predictors on
identical rows and identical folds, and compared by DeLong's paired estimator in the O(n log n)
midrank form of [9, 10]. The paired estimator is required rather than optional: the two score
vectors are fitted on overlapping training data and evaluated on the same rows, so they are
strongly correlated and treating them as independent understates the variance of the difference.
Confidence intervals on a single dataset's contribution are Wald intervals on the DeLong
standard error.

Both arms are evaluated out of fold. This is not a formality: an in-sample composition AUROC
compared against an out-of-fold model AUROC flatters composition enough to reverse the
comparison on individual datasets, reading as composition beating the model when it does not.

## Models

**4-mer logistic regression.** Raw 4-mer counts over the RNA sequence (98 per window, 256
columns) into L2-penalised logistic regression (C = 1). One model per fold, out-of-fold score
taken as the linear predictor. The k-mer order is 4; it is unrelated to the 5 in 5-fold
cross-validation, and the contrast is positive at every order from 3 to 6.

**Convolutional network.** A DeepBind-style architecture [11]: convolution (4 to 16 channels,
width 12), ReLU, max-pool 4, convolution (16 to 32, width 8), ReLU, global max-pool over
position, dense 32 to 64, ReLU, dropout 0.5, dense 64 to 1. 7,089 parameters, all trained from
scratch. Global max-pooling makes the model position-invariant, which is a requirement rather
than a convenience here: exploratory analysis placed the discriminative signal about 15 nt off
centre, varying by protein.

**SpliceBERT.** The pretrained nucleotide language model of [12], `multimolecule/splicebert`,
fully fine-tuned; 19.78 M parameters, all trainable. The classification head is a mean pool over
real nucleotide tokens, masking classification, separator and padding tokens, then dense to 128,
ReLU, dropout 0.3, dense to 1. Weights are baked into the container image at build time, so
"which weights" has the same answer as "which image".

Both neural models: AdamW, weight decay 0.01, batch size 32, binary cross-entropy, single
precision, at most 12 epochs with early stopping on validation AUROC at patience 4, best
checkpoint restored before test scoring. Learning rate 1e-3 for the head and 3e-5 for a
fine-tuned encoder. Nominal epochs overstate the training actually performed: the median
best epoch for SpliceBERT is 3 and most runs stopped by epoch 6 or 7, whereas the convolutional
network typically ran the full 12. Initial weights were drawn before the seed was set, so
initialisation is unseeded across the 940 committed fold-runs; the measured cost is a
per-dataset standard deviation of 0.006 to 0.010 in the nested contribution, inducing about
0.001 on a 94-dataset mean against a reported interval half-width of about 0.008.

## Inference

Because fifteen of the 79 proteins contribute two datasets each and the within-protein
correlation of the primary contrast is 0.92, every headline interval is a percentile interval
from a bootstrap that resamples **proteins** and takes all datasets of each sampled protein,
4,000 draws. Dataset-level intervals are 1.05 to 1.23 times narrower and are not reported as
headline values.

Per-dataset significance counts are additionally reported at a design effect of 1.35 applied to
the DeLong standard error. That factor is the product of two quantities measured on these data:
1.10 for DeLong conditioning on two score vectors as fixed when both are in fact fitted, and
1.23 for treating spatially clustered windows as independent, estimated by a 10 kb block
bootstrap. It moves the count of datasets where the 4-mer significantly helps from 80 of 94 to
72 of 94.

## Reproducibility and integrity

All published values are reproduced offline by a single command against committed evidence. The
harness runs 614 numeric assertions against a frozen expectations file, and the number of
assertions that ran is itself asserted, so a check cannot silently skip. Two assertions are
stronger than regression gates: 285 published AUROCs are rebuilt from committed per-example
scores to a maximum absolute difference of 2.2e-16, and the headline contrast is rebuilt from
raw sequence to 1.2e-06. Per-window out-of-fold scores for all three model classes are committed
(2,650 files), so every cell of the model-class comparison is recomputable from the repository
alone rather than being asserted against a hand-entered table.

Three defences against platform-dependent results are worth naming because they are otherwise
invisible: dinucleotide distances are computed on integer counts, candidate ties are broken by
genomic position, and the panel sort is stable. Each was added after a specific observed
discrepancy between CPU architectures or between runs.

Two limitations of the harness are stated in the Discussion: it verifies values rather than the
provenance of the code that produced them, and we found one case where that distinction
mattered.
