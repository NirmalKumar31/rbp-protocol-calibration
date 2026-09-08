# Changelog

## 1.0.0, the archived preprint release

Tagged `v1.0.0` and archived to Zenodo. This is the snapshot the manuscript cites and the one a
version DOI resolves to. The version identifies the snapshot; it makes no support promise, and
`SECURITY.md` still makes none.

The state the manuscript describes. 1156 numeric assertions pass offline against committed
tables; 882 collected tests; the paper builds from a clean export with zero LaTeX warnings, zero undefined references and zero over- or underfull boxes. "Warning-clean" used to be the phrasing and an audit objected, correctly, because a fresh log then carried one underfull \hbox: a typesetting diagnostic rather than a LaTeX Warning, but not nothing. That box was in methods.tex, where an unbreakable \texttt path forced TeX to stretch the line before it; the path now carries the same \allowbreak the manuscript already used in three other places, and both logs are at zero. Both counts are stated rather than one being folded into the other.

### The finding

Across 94 paired ENCODE eCLIP datasets, holding the model class, the peak set, the
chromosome-to-fold map and the estimator fixed and varying how negative windows are built moves a model's
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

### What the eleventh review round changed

An external audit of the frozen candidate found no new result-changing analytical defect and
eight release defects. All eight are closed:

- **The Discussion quoted a withdrawn result as the current one.** It gave the
  chromosome-blocked cross-fitted external cell as 1.68 (1.41 to 2.02), the value the shared-map
  repair replaced with 1.67 (1.40 to 2.00), and in the same sentence gave the cross-fold
  near-neighbour fraction as 0.0036, the strand-blind figure the strand repair withdrew. The
  abstract, the Results table and the Results prose were all correct; only the Discussion was
  not. Nothing caught it, because both old values remain legitimately present in
  `docs/EXTERNAL_CORRECTION_1.md`, which preserves them on purpose, and both are stated to two
  decimals, which `audit_manuscript.py` deliberately does not trace.
- **All four external cells are now gated on every surface that states them** — abstract,
  Results table, Results prose, Discussion and README — against `external_sensitivity.csv`,
  point and both bounds, recomputed rather than pinned.
- **The README's assertion split was stale and did not add up.** It said 1015 of 1156 belonged
  to this paper where the verifier says 1018, and 1018 + 136 is 1154. The missing pair is the
  harness checking itself; all three parts are now derived from `verify_summary.csv` and gated.
- **The documented clean install could not run the checks it then invokes.** `pip install -e .`
  omits pytest, Ruff and pypdf, which are in the `dev` extra, and a constraints file pins what
  is being installed rather than adding anything. The README now installs `.[dev]` and
  `.[dev,neural]`.
- **Two preflight steps could report success without checking anything.** The shell-syntax step
  ran `git ls-files` and, outside a checkout, silently looped over zero files: in a `git archive`
  extraction holding every tracked shell script it printed OK having parsed none of them. And the no-pdflatex
  branch printed the word SKIPPED with a bare `printf`, so the flag stayed clear and the run
  still ended on PREFLIGHT CLEAN with the PDF gate unrun.
- **`run.sh` and `docs/REPRODUCE.md` called the default path the whole study.** It is the core
  pipeline: the dinucleotide arm end to end, with the GC and bias-aware neural sweeps consumed
  from committed per-window scores. The detailed caveats were already correct; the headline was
  not.
- **Terraform and `config/params.yaml` named different GCP projects.** The variable default and
  `terraform.tfvars.example` both carried the previous study's project, so accepting the
  defaults would have put the state in one project and every resource in another. The default is
  removed, the example matches, and `run.sh` refuses to init on a mismatch.
- **GitHub reported the repository licence as NOASSERTION**, because the data and results
  licensing was appended to the MIT text. `LICENSE` is now the unmodified MIT licence and that
  section moved verbatim to `NOTICE`, with the two CC BY files repointed and the cross-reference
  test widened to follow them.

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
