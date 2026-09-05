# Changelog

## 0.9.0 — release candidate for the preprint

Not tagged. `1.0.0` is reserved for the commit that is archived and given a DOI.

The state the manuscript describes. 998 numeric assertions pass offline against committed
tables; 727 collected tests; the paper builds warning-clean from a clean export.

### The finding

Across 94 paired ENCODE eCLIP datasets, holding the model, the peak set, the chromosome-blocked
folds and the estimator fixed and varying only how negative windows are built moves a model's
measured nested contribution 5.42-fold, while its apparent AUROC moves the opposite way.

### What the last review round changed

An external audit of the release found 60 items. The substantive ones:

- **The estimator's floor is the outer-fold information route, not conditioning.** Rerunning
  with the route closed cuts it by at least 95% and recovers the zero that theory requires to
  within 5e-4. The Results had attributed it to conditioning; corrected, with the measurement
  in a new section and a gate on it.
- **The bias is not one-directional.** The Methods claimed it "can only help the score column".
  For the 4-mer, closing the route *raises* the contribution in all three arms. Withdrawn.
- **The bias-aware arms had no committed panel.** Their membership came from directory-listing
  an uncommitted store. Now committed and cross-checked against the per-window evidence.
- Four overclaims narrowed to what the design supports: "invariant to fold size", "could only
  widen", "unbiased by construction", "independent benchmark".

### Infrastructure

- Numeric audit widened to `paper.tex` — the title and abstract had never been checked — and to
  the five release documents.
- `scripts/release_consistency.py`: every count a document states about the release is derived
  from the release and gated.
- CI now runs the verifier, the release check, ruff, shell syntax and the LaTeX build. It ran
  none of them before.
- Ruff 193 violations to 0.
- Modal cost guard fixed: it matched one app name while three of four sweeps ran under others,
  counted one arm of four, and failed open.
- `CITATION.cff`, complete package metadata, per-directory data licences, a column dictionary,
  one authoritative cost table.

### Known limits, stated rather than closed

- Cross-fitting is measured for the k-mer classes only; the CNN and SpliceBERT would need four
  times the GPU sweep.
- One negative draw, one fold partition, unseeded neural initialisation. Variability is
  quantified in the Limitations and is not propagated into the headline intervals.
- Sequence homology across folds is audited by exact 32-mer sharing, not by identity clustering.
- The abstract is over length for most venues.
