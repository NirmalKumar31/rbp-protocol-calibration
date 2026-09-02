# Introduction

Predicting where an RNA-binding protein (RBP) binds from sequence alone is a standard supervised
problem, and the field's progress is measured almost entirely by AUROC on held-out windows.
Positives come from a crosslinking assay, usually eCLIP [1, 2], and the negatives are
constructed. Horlacher and colleagues [3] made the consequences of that construction explicit:
across eleven RBP prediction methods and hundreds of experiments, the apparent performance of a
method depends on how its negative set was built, and they proposed a bias-aware alternative in
which negatives are drawn from other RBPs' binding sites rather than from genomic background.
That the number moves is therefore established. What has not been established is what the number
is a function of, how far it moves once a composition baseline is accounted for, and whether any
reported quantity survives the change.

The question matters because a raw AUROC conflates two things a reader cannot separate: how much
signal the model found, and how easy the discrimination was made. A window set whose negatives
are sampled uniformly from the genome differs from its positives in base composition alone,
and a model that detects "G-rich, C-poor" will score well without having learned anything
sequence-specific. The standard mitigation is to match negatives on composition, most often GC
content. But GC content is one linear functional of the sixteen-cell dinucleotide simplex, so
matching on it leaves fifteen degrees of freedom free, and the residual compositional difference
is available to any model. This is the same structure that Grimm and colleagues [4] identified in
variant-effect benchmarks, where the way negatives and positives are assembled produces
circularity that inflates apparent method performance and reorders methods, and it is a
recognised failure mode of genomic machine-learning evaluation more generally [5, 6].

The natural response is to report not the model's AUROC but what it adds over an explicit
composition baseline. That quantity, the incremental or nested contribution, has a long history
in risk prediction [7, 8], and it has an exact meaning: fit the composition features alone and
then with the model's score, on the same rows and the same folds, and take the difference in
out-of-fold AUROC. It is the quantity a reader actually wants when asking what a model
contributes, and it is the quantity we measure here.

Our contribution is a calibration. We take 94 ENCODE eCLIP datasets, hold the model, the
positives, the chromosome-blocked folds and the estimator fixed, and build negatives three
ways: matched on GC content, matched on all sixteen dinucleotide frequencies, and drawn from
other RBPs' binding sites in the same cell line following the bias-aware recipe of [3]. Under
each protocol we refit the composition baseline on that protocol's own windows, so the baseline
is never imported across designs. Four results follow.

First, apparent difficulty and measured contribution move in **opposite** directions.
Dinucleotide matching lowers apparent AUROC in every one of the 94 datasets while roughly
tripling the measured contribution, and the bias-aware protocol, which the field reads as the
most principled of the three, turns out to be the **easiest** discrimination and yields the
**least** measured contribution. A reader who infers from a lower score that a benchmark is more
demanding, or from a higher score that a model is better, can be wrong in a direction that no
amount of care with the model can fix.

Second, the range is large in the units that matter. The same 4-mer on the same positives
measures 0.066, 0.027 or 0.012 AUROC of contribution depending only on the negatives, a
5.4-fold span, and the absolute spread of 0.054 is about three times the 0.019 that separates a
4-mer logistic regression from a convolutional network on our own panel. We are deliberate about
which comparison this licenses. Protocol choice exceeds the difference between adjacent
conventional model classes, so a benchmark that changes its negatives can reorder such methods.
It does not exceed every architectural difference: the step from the convolutional network to a
fine-tuned nucleotide language model is 0.103 on the same panel, about twice the protocol
spread. The claim is that protocol effects are of the same order as, and often larger than, the
effects these benchmarks are built to resolve, not that they dominate all of them.

Third, the composition baseline each protocol leaves behind accounts for most of the variation,
and we report that mechanism rather than defending against it. Given the baseline, knowing which
protocol produced it adds about one percent of variance; given the protocol, knowing the
baseline adds eleven. The relationship is not uniform, however: it is a property of
composition-matched negatives and is absent for the bias-aware arm, so there are two protocol
families rather than three interchangeable protocols, and a residual difference between families
survives at matched baseline. That family structure replicates on an independent benchmark.

Fourth, no rescaling supplies a transportable fix. AUROC is compressive near one, so part of the
span is arithmetic, and eight standard monotone reparameterisations reduce it. None removes it:
the smallest span over the eight is 2.00-fold, far outside the range expected if the three
protocols had equal true means. An exponent that does equalise our three protocols exists, but
it is fitted, and the exponent fitted to an independent benchmark differs from ours by a factor
of 2.4, leaving that benchmark essentially where it started. The absence of a protocol-free
measure is thus a statement about transportability, which is testable, rather than about a floor,
which is not.

We therefore recommend a two-number report: state the composition-only AUROC obtained under the
same protocol beside every headline AUROC. We are explicit about what this does and does not
buy. It makes the problem visible and it is the coordinate least sensitive to protocol among
those we examined, but on an independent benchmark it does not improve cross-protocol agreement,
so we do not claim it makes contributions comparable. It does not; we recommend not comparing
them.

Throughout, we use "measured contribution" and never "true contribution". The three protocols do
not estimate a common quantity, so words like hidden, concealed, inflated or proper would
presuppose a value this work argues does not travel between benchmarks. The register is
calibration.
