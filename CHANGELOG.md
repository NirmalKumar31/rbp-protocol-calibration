# Changelog

## 1.0.0, preprint release

This version identifies the snapshot prepared for archival and makes no support promise;
`SECURITY.md` states the scope.

### The finding

Across 94 paired ENCODE eCLIP datasets, holding the model class, the peak set, the
chromosome-to-fold map and the estimator fixed and varying only how negative windows are built,
a model's measured nested contribution moves 5.42-fold under the conventional two-stage
estimator and 4.84-fold under the cross-fitted estimator that is this paper's primary. Its
apparent AUROC moves the opposite way: the protocol that looks hardest yields the largest
contribution.

The ordering holds across three model classes, five estimands, four class ratios and every
baseline order tested, and replicates on 135 held-out datasets from an independently constructed
benchmark. The inverse direction does not replicate there, and the paper says so.

### What is in the release

- The manuscript and supplement, and the LaTeX that builds them.
- Every result table the paper cites, with a column dictionary and a provenance record
  classifying each table as raw-reproducible, evidence-recomputable or frozen.
- Per-window out-of-fold scores for all three model classes, so the model-class comparison is
  recomputable from the repository alone rather than asserted against a summary.
- The analysis code, the cloud definitions that ran it, and the container images.
- 1156 numeric assertions that check the published values offline against committed evidence,
  in one command and with no cloud account.
- Four pre-specified protocols under `docs/`, cited by name in the paper, so the claims that
  rest on them can be read rather than taken on trust.

### What the verification does and does not establish

The suite is a regression and provenance gate. Most assertions compare a reported value with
committed evidence from the same run. Two are stronger: 285 AUROCs are rebuilt from per-window
scores, and the headline contrast is rebuilt from raw sequence. The manuscript-number trace can
identify a value no table supports, but not whether a supported value is attached to the right
claim. Scientific validity rests on the design and the limitations, not on the gate count.

### Known limits, stated rather than closed

- Cross-fitting is measured for the k-mer classes only. Closing the same route for the
  convolutional network and SpliceBERT would cost about $115 in GPU time and was not run, so
  their spans are exploratory two-stage estimates.
- One negative draw, one fold partition, unseeded neural initialisation. Each is quantified in
  the Limitations and none is propagated into the headline intervals.
- Sequence homology across folds is audited by exact 32-mer sharing, not by identity clustering.
- The panel is ENCODE eCLIP in two cell lines. The calibration is for eCLIP-derived panels.
