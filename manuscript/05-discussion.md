# Discussion

The practical consequence of these measurements is narrow and actionable. When reporting how
well a sequence model identifies RBP binding sites, report the composition-only AUROC obtained
under the same protocol beside the headline AUROC, and do not compare contributions measured
under different negative-set protocols. The first is cheap: it is one extra logistic regression
on nineteen features. The second is a restriction on inference rather than an addition to it,
and it is the part that costs something, because much of the field's cross-paper comparison is
of exactly that kind.

Three things follow that are worth separating from the recommendation itself.

**Apparent difficulty is not a proxy for benchmark quality.** The intuition that a lower score
on a harder negative set indicates a more honest evaluation is what the reversal breaks. Moving
to dinucleotide-matched negatives lowers apparent AUROC in every dataset we examined while
raising the measured contribution nearly threefold, and the bias-aware protocol, which removes
real artefacts and is well motivated on those grounds, produces the *easiest* discrimination and
the *smallest* measured contribution of the three. Nothing about the bias-aware design is wrong;
what is wrong is reading its lower or higher scores as information about how much a model has
learned. Different proteins bind compositionally different sites, so separating one protein's
sites from another's is largely a composition task, and a protocol that makes that the task
leaves a sequence model little to add.

**There are two protocol families, not a continuum of stringency.** The relation between the
composition baseline a protocol leaves and the contribution a model then shows is strong for
composition-matched negatives, in our data and independently in [4]'s, and absent for negatives
drawn from other proteins' sites. A residual difference between the families survives matching
on the baseline. So the useful taxonomy for a benchmark builder is not "how strict is the
matching" but "what kind of object is the negative", and results are comparable within a family
in a way they are not across families.

**A quantity's protocol dependence can be reduced by fitting and not by choosing.** Eight
standard reparameterisations shrink the span without closing it, and the best of them, dividing
by the baseline's headroom, is the two-number report expressed as one number. But an exponent
that closes the span entirely exists, and it does not transport: the exponent fitted to our
three protocols leaves an independent benchmark almost exactly where it started, and the
exponent fitted to that benchmark is 2.4 times ours. This is the sense in which there is no
protocol-free measure of what a sequence model contributes. It is not that no transformation
works; it is that the transformation that works has to be refitted per benchmark, which is
precisely what a comparable measure cannot require.

**On what "adds over composition" means, and it is a real qualification.** The composition
baseline stops at order two, which is what a dinucleotide-preserving shuffle would have
preserved and is a defensible place to stop. It is nonetheless a choice, and it does substantial
work. Refitting the baseline at order three, adding the 63 trinucleotide frequencies, removes
about four fifths of the 4-mer's measured contribution and 87% of the headline contrast, and on
the bias-aware arm the median dataset retains nothing at all. **For a bag-of-k-mers model, most
of what we and others report as "what the model adds over composition" is one further order of
composition, not motif recognition and nothing positional.** Any paper reporting an increment
over a composition baseline is reporting a quantity indexed by that baseline's order, and almost
none state which order they used.

Two things temper that, and one of them is a result in its own right. The protocol dependence is
not an artefact of where the baseline stops: at an order-three baseline the contrast remains
+0.0067, +0.0309 and +0.0545 for the three model classes, all clear of zero, and the span across
protocols does not collapse. And the collapse is specific to models whose features are nested
inside the baseline's span. Over a trinucleotide baseline the 4-mer's contribution is positive in
only 65 of 94 datasets while SpliceBERT's is positive in 94 of 94, and once the raised baseline's
reduced headroom is accounted for the 4-mer loses about 1.4 times what SpliceBERT does. A
benchmark that wants to separate model classes rather than rank them by how much short-range
composition they absorb could therefore raise the baseline's order. We state the circularity that
recommendation carries: trinucleotide counts are linear aggregates of 4-mer counts, so an
order-three baseline is structurally unfavourable to k-mer models, and it separates models by how
much of their signal lies outside short-range composition rather than by capacity.

**Relation to prior work.** [4] established that negative-set construction changes apparent
performance across many RBP methods and proposed the bias-aware alternative. This work does not
contest that and depends on it: their released negative sets are our only external validation.
What we add is the calibration. The quantity that moves is not just apparent AUROC but the
increment over an explicit baseline; it moves in the opposite direction to apparent difficulty;
it moves by a factor of five, comparable to or larger than the between-method differences these
benchmarks exist to resolve; the composition baseline rather than the protocol label is what it
tracks; and no reparameterisation makes it comparable across benchmarks. The same structural
problem is well documented in variant-effect prediction [6], where benchmark assembly produces
circularity that both inflates and reorders methods, and the general point that
genomic machine-learning evaluations are sensitive to how the negative class is defined has been
made before [6, 7]. Our contribution is to quantify it for one well-defined estimand in one
well-covered assay, and to show which part of it is arithmetic and which is not.

## Limitations

**The sign of the primary contrast is design-implied.** GC matching constrains one of the
composition baseline's fifteen degrees of freedom and dinucleotide matching constrains all
fifteen, so the direction of that contrast follows from the design and only its magnitude is
informative. The bias-aware arm is not subject to this argument, since it constrains none of the
fifteen and still has the highest baseline.

**Every magnitude is indexed to a mono- and dinucleotide baseline.** See the Discussion above.
The protocol dependence is not so indexed; the magnitudes are.

**Achieved matching is not exact matching.** The GC matcher's operative acceptance bound is 0.15
rather than its nominal 0.05, holding the nominal tolerance for 94.8% of pairs, and the
dinucleotide matcher reduces composition distance 2.27-fold rather than to zero. Degrees-of-
freedom statements are exact for the design and approximate for the data.

**Three protocols, two cell lines, three model classes, one assay.** The three-protocol results
use the 4-mer only, because the neural models were never trained on the bias-aware arm; the
multi-model results use two arms only. The two halves of the argument therefore do not meet in a
single cell, and extending the neural sweep to the third protocol is the most valuable single
addition to this work.

**Provenance limitations on the neural arms.** For 20 of 94 datasets in the dinucleotide arm the
committed neural scores came from a partition that is not chromosome-grouped; we report the 74
correctly partitioned datasets as primary and the full panel as a sensitivity, with shifts inside
the interval half-widths. Initialisation was unseeded across the 940 fold-runs, at a measured
cost of about 0.001 on a panel mean. The convolutional rung was trained on a different accelerator
from the transformer rung. And the neural scores entered the nested fit on the probability scale
while the 4-mer's entered as a log-odds, which understates the neural contributions by up to
0.0033 and moves the contrasts by at most 0.0004.

**Positive sets differ slightly between arms.** Window construction is identical, but a positive
is dropped when its matcher finds no acceptable negative, and the two matchers fail on different
windows. Measured across the panel, the two arms' positive sets have a Jaccard similarity with
median 0.9972 and minimum 0.9237, and are exactly identical in 10 of 94 datasets
(Supplementary Table S9). The negatives are therefore almost, and not exactly, the only
difference between arms.

**No significance filter on positives and no expression filter on negatives.** Every called peak
is a positive, so the positive set includes weak peaks that a study applying ENCODE's standard
thresholds would drop. The expression question is bounded by the control reported above rather
than removed. Region labels are 6.1% to 8.7% wrong on the negative side by our own audit, which
affects the region-matching step asymmetrically.

**Two limitations of the verification harness itself**, stated because the harness is otherwise
easy to over-read. It asserts values, not the provenance of the code that produced them, and the
fold-partition defect above is exactly the class of error that distinction permits: it was found
by inspection, not by any of the 614 assertions. And its manuscript audit only checks numbers
written to three or more decimal places, so a percentage written to one decimal is invisible to
it and was checked by hand.

**Inference is asymptotic in places.** Confidence intervals on a single dataset's contribution
are Wald intervals on the DeLong standard error rather than bootstrap intervals, and the
published Firth coefficient intervals resample rows rather than genes, which is narrower than
gene clustering would give.
