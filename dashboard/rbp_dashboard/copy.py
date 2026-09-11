"""Prose the dashboard shows: term definitions, and the caveats that must travel with figures.

Kept here rather than inline so that a claim appears once. The wording follows the manuscript;
where it compresses, it compresses toward the weaker statement.
"""

GLOSSARY = {
    "AUROC": (
        "Area under the ROC curve. The chance that a randomly chosen bound window scores above "
        "a randomly chosen unbound one. 0.5 is coin-flipping, 1.0 is perfect."
    ),
    "negative set": (
        "The windows called 'not bound'. They are chosen, not observed, and this study is about "
        "how much that choice changes the answer."
    ),
    "composition baseline": (
        "A model given only 19 numbers describing nucleotide composition: how much A, C, G and T, "
        "GC fraction, and dinucleotide frequencies. No sequence order, no learning about motifs."
    ),
    "nested contribution": (
        "AUROC of (composition features + the model's score) minus AUROC of composition alone. "
        "What the model adds once plain composition is already accounted for."
    ),
    "apparent AUROC": (
        "The score a paper would print as its headline: composition features plus the model's "
        "score, with no baseline subtracted."
    ),
    "GC-matched": (
        "Unbound windows picked so their overall G+C fraction resembles the bound windows'."
    ),
    "dinucleotide-matched": (
        "A stricter match: pairs of adjacent bases are matched too, not just single-base content."
    ),
    "bias-aware": (
        "Unbound windows drawn to control the technical biases of the assay itself, rather than "
        "matched on composition."
    ),
    "two-stage estimator": (
        "Fit the model, then measure what it adds. The literature surveyed here does this. It "
        "lets information from the evaluation folds leak into the score being evaluated."
    ),
    "cross-fitted estimator": (
        "The model's score for each window is produced by a model that never saw that window's "
        "fold. This closes the leak. It is the primary estimator in this paper."
    ),
    "estimator floor": (
        "What an estimator reports when the true answer is known to be zero. Anything above zero "
        "is the estimator measuring itself."
    ),
    "protein-clustered bootstrap": (
        "Resampling whole proteins rather than single datasets, so two datasets for the same "
        "protein cannot be treated as independent evidence."
    ),
    "span": (
        "Largest divided by smallest across the three negative-set protocols. A span of 1.0 "
        "would mean the protocol made no difference."
    ),
    "held-out benchmark": (
        "135 datasets from a separate published deposit, sharing no dataset with the 94 analysed "
        "here, used to test whether the finding survives outside this pipeline."
    ),
}

# Every figure comparing model classes must carry this. The 4-mer result is cross-fitted and
# primary; the neural results are two-stage and exploratory, and conflating them would overstate
# what the study establishes.
ESTIMATOR_CAVEAT = (
    "The 4-mer result is reported under both estimators, and its cross-fitted value is the "
    "paper's primary estimand. The CNN and SpliceBERT spans are computed with the conventional "
    "two-stage estimator only and are exploratory. They were not cross-fitted, so they carry the "
    "same outer-fold channel the 2-mer test below shows to be non-zero."
)

PANEL_CAVEAT = (
    "95 datasets were selected for the panel. 94 carry all three negative-set protocols and are "
    "the denominator for every three-protocol comparison here."
)

NOT_ONLY_NEGATIVES = (
    "The negative-set construction is what varies, but it is not literally the only thing that "
    "changes. The model and the composition baseline are refitted inside each protocol, and each "
    "protocol keeps the positives its matcher could pair, so retained positive subsets differ "
    "slightly between arms. The manuscript quantifies that difference and shows it does not "
    "account for the result."
)

REPRO_CAVEAT = (
    "This is verification, not complete raw-to-result reproduction. The release re-derives each "
    "committed table from cached intermediates and checks it against the frozen copy. Regenerating "
    "everything from raw ENCODE reads is a separate and longer path, and no claim of bit-for-bit "
    "end-to-end reproducibility is made."
)

NOT_MONITORING = (
    "A static reader over committed CSV files. It is not production monitoring, it is not "
    "connected to any cloud service, and it does not run analyses or recompute estimates."
)
