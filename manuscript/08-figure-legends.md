# Figure and table legends

Every legend states n, the resampling unit and what any interval covers, because those three
things are what a reader cannot recover from the plot. Panel letters match the in-axes labels.
Source table is given for each figure so a reader can go from picture to number to evidence.

## Main figures

**Figure 1. The three protocols, and the reversal.** *(f10_three_protocols)* Same 94 ENCODE
eCLIP datasets, same positives, same 4-mer model, same chromosome-blocked folds and same
estimator throughout; only the construction of the negative windows differs. **(a)** Composition
-only AUROC, apparent AUROC and nested contribution under each protocol. **(b)** Nested
contribution per dataset, dinucleotide-matched against GC-matched, paired; points above the
diagonal are datasets where the model contributes more under dinucleotide matching (88 of 94).
**(c)** The same comparison for the bias-aware arm against GC-matched, where the contribution is
lower in 64 of 94. Note that apparent AUROC and contribution move in opposite directions between
panels: dinucleotide matching lowers apparent AUROC in 94 of 94 datasets while raising the
contribution. Source: `three_arm_contrast.csv`, `three_arm_per_dataset.csv`.

**Figure 2. No transportable rescaling.** *(f11_scale_sweep)* Fold range in measured contribution
across the three protocols, under eight monotone reparameterisations of the same 282
protocol-dataset cells. Points are the range; bars are 95% percentile intervals from a bootstrap
over datasets, 2,000 draws. The dashed line at 1.0 marks protocol independence. Somers' D is an
affine map of AUROC and returns the raw value by construction, serving as a control on the sweep.
Normalising by the excess over chance inflates the range to 18.2 and is shown because it is the
transformation a reader is most likely to invent. The minimum over the eight, 2.00 [1.67, 2.46],
is attained by dividing by the baseline's own headroom. Two qualifications are in the text: the
comparison point is not 1.0 but the equal-means null (median 1.06, 95th percentile 1.18), and an
exponent that does equalise the three protocols exists but does not transport. Source:
`scale_sweep.csv`, `transport_check.csv`.

**Figure 3. It is the composition baseline, not the protocol label.** *(f12_protocol_or_baseline)*
**(a)** All 282 protocol-dataset cells: nested contribution against the composition-only AUROC
the protocol left behind, coloured by protocol. Pooled Spearman −0.600 (p = 6e-29). **(b)** The
natural experiment. Datasets are split by whether the bias-aware protocol raised or lowered the
composition baseline relative to GC matching, and the difference in contribution is shown for
each stratum. Points are stratum means; bars are 95% percentile intervals from a bootstrap over
datasets, 4,000 draws; n is 67 and 27. The deficit reverses where the baseline falls, so the sign
follows the baseline rather than the protocol label. Source: `three_arm_per_dataset.csv`,
`protocol_or_baseline.csv`.

**Figure 4. The effect is not a property of the model class.** *(f9_deep_contrast)* Nested
contribution in each arm for a 4-mer logistic regression, a 7,089-parameter convolutional network
and a fully fine-tuned 19.78 M-parameter SpliceBERT, on identical rows and folds, n = 94.
Intervals are 95% percentile intervals from a bootstrap over the 79 **proteins**, taking all
datasets of each sampled protein, 4,000 draws. The contrast is present for all three classes; on
the ratio scale the ordering reverses, so the contrast does not grow with capacity. For the two
neural models the primary panel in the text is the 74 datasets whose committed scores are
chromosome-partitioned; see Limitations. Source: `deep_contrast.csv`,
`deep_contrast_per_dataset.csv`, `fold_integrity.csv`.

**Figure 5. The calibration replicates on an independent benchmark.** *(f14_external_validation)*
Horlacher *et al.*'s released negative sets, their peak calling and their fold assignments, over
the 45 datasets shared with our panel; only the measurement is ours. **(a)** Nested contribution
under their negative-2 (other RBPs' sites) against their negative-1 (transcript background),
paired per dataset; 36 of 45 fall below the diagonal, a 2.38-fold ratio of panel means. **(b)**
The mechanism, which is what actually travels: within-arm Spearman correlation between the
composition baseline and the contribution, for the two composition-matched arms in our data, the
composition-matched arm in theirs, and the other-RBPs'-sites arms in both. The gradient is present
in the composition-matched family in both benchmarks (−0.54, −0.46, −0.64) and absent in the
other-RBPs'-sites family in both (−0.12 and −0.19, neither significant). Source:
`horlacher_per_dataset.csv`, `three_arm_per_dataset.csv`.

**Figure 6. The recommendation, and where it stops working.** *(f15_recommendation)* **(a)** On
our data, cross-protocol rank agreement between the contributions a reader would compare, before
and after dividing by the baseline's headroom, for each of the three protocol pairs. Agreement
improves in 3 of 3 pairs, but only the GC-versus-dinucleotide improvement has a protein-clustered
interval clear of zero (+0.051 [+0.007, +0.115]); the others are marked n.s. **(b)** The same test
on the independent benchmark of Figure 5, n = 45. Rank agreement falls and disagreement rises
under normalisation. Those are the two criteria pre-registered as falsifying the recommendation,
and both fire; the change in rank agreement is −0.050 [−0.222, +0.140], so this is a failure to
replicate rather than a refutation. Panels (a) and (b) are shown together deliberately: a figure
presenting only (a) would be the strongest available misrepresentation of the paper's own
evidence. Source: `recommendation_works.csv`, `transport_check.csv`.

## Main table

**Table 1. The model ladder and the panel.** Out-of-fold AUROC on identical chromosome folds,
n = 95 datasets: composition alone 0.628, 4-mer 0.688, convolutional network 0.706, fine-tuned
SpliceBERT 0.809. Included as a methods table rather than a result, both to establish that the
models behave as expected on this panel and to give the scale against which the protocol effects
in Figures 1 and 2 should be read: the gap between the 4-mer and the convolutional network is
0.019, against a protocol-induced spread of 0.054. Source: `matched_four_models.csv`.

## Supplementary figures

**Figure S1. Panel composition.** *(f0_panel_overview)* **(a)** Study panel size distribution
against the candidate pool it was drawn from; the panel spans the pool's 0th to 99th percentile,
which matters because AUROC correlates with dataset size at r = 0.53 to 0.67 here. **(b)** Split
by cell line and count of proteins assayed in both. **(c)** The composition-only AUROC each
protocol leaves, as three distributions; their near-disjointness is why protocol and baseline
cannot be separated by a pooled analysis. Source: `panel_summary.csv`, `candidate_sizes.csv`,
`three_arm_per_dataset.csv`.

**Figure S2. The two-arm contrast, paired.** *(f1_cost_of_matching)* n = 94 paired datasets.
Source: `cost_of_matching.csv`.

**Figure S3. The model ladder.** *(f2_four_models)* n = 95, identical folds. Source:
`matched_four_models.csv`.

**Figure S4. The strand control against its matched placebo.** *(f3_strand_placebo)* n = 40
datasets over 40 distinct proteins, so no clustering correction is needed. Twenty placebo seeds;
intervals are 95% percentile intervals over datasets. Source: `strand_placebo.csv`.

**Figure S5. Cross-cell-line replication.** *(f4_replication)* The 15 proteins assayed in both
K562 and HepG2. Source: `multidonor_pairs.csv`, `robustness.csv`.

**Figure S6. Robustness to k.** *(f6_k_sweep)* The contrast is positive at every k from 3 to 6,
and positive at every k individually in 82 of 94 datasets. Source: `k_sweep.csv`.

**Figure S7. Where the composition baseline stops.** *(f13_baseline_order_models)* Share of each
model's contribution surviving an order-3 (mono + di + tri) composition baseline, both arms,
n = 94. Intervals are protein-clustered, 2,000 draws. The shares differ mainly because the totals
differ: the absolute amount absorbed is near-constant across model classes, and after correcting
for the headroom the raised baseline removes, the 4-mer loses about 1.4 times what SpliceBERT
does. See the Discussion. Source: `baseline_order_models.csv`.

**Figure S8. Effect modification by dataset size.** *(f5_size_modification)* Spearman 0.307
(p = 0.0026) between dataset size and the reported contrast. Source: `cost_of_matching.csv`.

## Supplementary tables

**Table S1. The panel, joined to ENCODE.** One row per dataset: protein, cell line, ENCFF file
accession, ENCSR experiment accession, replicate count, pair count, whether the dataset carries
all three protocols, and each protocol's composition baseline and nested contribution. 95
datasets, 79 proteins, 95 distinct experiments. 94 datasets carry all three protocols and are the
unit of every three-protocol result; one is ladder-only. Source: `supplementary_table_s1.csv`.

**Table S2. Achieved match quality.** Per arm and per dataset: median, p90 and maximum absolute
GC gap and dinucleotide L1 distance between each positive and its assigned negative, over all
456,734 pairs per arm. Source: `match_quality.csv`, `match_quality_per_dataset.csv`.

**Table S3. Fold integrity of the committed scores.** Per dataset, arm and model: agreement
between the score file's fold assignment and the study folds, the maximum number of chromosomes
in any one fold, and the fraction of rows with a same-strand genomic neighbour within 1 kb
assigned to a different fold. Source: `fold_integrity.csv`,
`fold_integrity_per_dataset.csv`.

**Tables S4 to S8.** Per-dataset tables underlying each supplementary section: the strand and
expression controls, the specification grid for the non-identified decomposition, the clustered
intervals, and the transform sweep. Sources: `strand_placebo_per_dataset.csv`,
`expression_control_per_dataset.csv`, `protocol_identification.csv`, `cluster_intervals.csv`,
`scale_sweep.csv`.

---

## Note on figure preparation

All figures are generated by `scripts/figures.py` from committed tables only; nothing is computed
in the plotting code that is not already in a table, so a figure cannot disagree with the text.
PDFs embed TrueType (Type 42) fonts rather than matplotlib's default Type 3, which this and
comparable journals reject. Rasters are 400 dpi. Colour is assigned by model or protocol identity
and never by rank, so a figure that drops a series does not repaint the others.

**Table S9. Positive-set overlap between arms.** Per dataset, the number of positive windows in
each composition-matched arm and the Jaccard similarity between them. A positive is dropped when
its matcher finds no acceptable negative, and the two matchers fail on different windows, so the
arms' positive sets are near-identical but not identical: median Jaccard 0.9972, minimum 0.9237,
exactly identical in 10 of 94 datasets. Source: `positive_set_overlap.csv`.
