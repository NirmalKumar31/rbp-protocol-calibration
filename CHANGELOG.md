# Changelog

## 1.0.4, the archive points at the paper

Metadata only. No code, no results, no manuscript text changes.

The Zenodo record carried one related identifier, a link to the GitHub tree, and nothing
saying which paper the deposit belongs to. `CITATION.cff` gained the preprint DOI in 1.0.3,
but Zenodo has no field that maps a CFF `preferred-citation` onto a related identifier, so
the archive stayed silent about the paper it supports.

- **`.zenodo.json` added**, declaring the posted preprint
  `10.20944/preprints202609.0883.v1` as `isSupplementTo`, resource type
  `publication-preprint`. Both vocabulary values were checked against Zenodo's live
  vocabularies rather than assumed.
- **It reproduces every field the record already had** (title, abstract, keywords, creator,
  ORCID), because Zenodo ignores `CITATION.cff` entirely once this file is present. Each
  field was compared against the published 1.0.3 record before release.
- **`version` and `license` are deliberately omitted.** Version comes from the git tag and
  licence from the LICENSE file, which is the only place the per-subtree licensing is
  recorded. Naming either here would create a copy nothing bumps.
- **A test pins the two files together.**
  `test_zenodo_json_and_citation_cff_cannot_drift` fails if the title, abstract, keywords,
  creator or ORCID diverge, if the preprint DOI leaves `.zenodo.json`, or if `version` or
  `license` appear in it. The duplication is the same defect shape the rest of that module
  exists to catch, so it is gated rather than trusted.

## 1.0.3, the posted preprint

No result changes. This version makes the repository, the Zenodo archive and the posted
paper say the same thing, which they did not between 1.0.2 and the Preprints.org posting.

- **The paper is posted.** Preprints.org, 2026-09-11, DOI `10.20944/preprints202609.0883.v1`,
  CC BY 4.0. `CITATION.cff` records it in `preferred-citation`, which previously carried no
  DOI at all and so left no machine-readable trace of where the work appeared.
- **Use of AI tools.** Methods gains a subsection naming the tools, what they were used for
  and who is responsible. Preprints.org requires the declaration; it belongs in the paper
  regardless, and the repository copy lacked it while the submitted copy had it.
- **Conflicts of Interest.** The declaration was titled "Competing interests" and sat before
  Ethics. It is now titled and positioned as the submitted version has it.
- **A false comment removed.** `CITATION.cff` asserted that the software title and the paper
  title "carry the same title on purpose". They differ, deliberately, and have since this
  record became `type: software`. The comment described a decision that had been reversed:
  a pointer that read confidently and was wrong, which is the defect class this paper is
  about.

`manuscript/paper.tex` and `manuscript/sections/methods.tex` are now byte-identical to the
source submitted to Preprints.org. Every published number is unchanged.

## 1.0.2, submission snapshot

No result changes. This version exists so that the PDF submitted to a preprint server is the
PDF that is archived, which was not true of 1.0.1.

- **Affiliation.** The title page named only "Northeastern University". It now names the
  department and college.
- **Funding statement.** It read "This research received no specific funding. I paid the
  cloud-computing costs." The second sentence is true and, beside a single author with a bare
  university line, reads as a personal project rather than as research. It is now the standard
  form. The full cost breakdown stays in `docs/COST.md`, where it is a strength.
- **Venue.** bioRxiv declined the submission at screening, without review, as a student
  research project. Comments and the README no longer name a target venue or claim a preprint.
- **README.** Rewritten for a reader rather than for an adversarial reviewer: 2,035 words to
  1,647, a first-person section explaining why the study pivoted from ClinVar variant scoring
  to negative-set calibration, the offline check moved near the top, and the findings table cut
  from a 114-word maximum cell to 35.

Every published number is unchanged. 1156 assertions pass, CI is green on all five jobs, and
`results/tables` is byte-identical to 1.0.0.

## 1.0.1, the DOIs

The 1.0.0 archive necessarily contained a paper with no DOI printed in it, because Zenodo
cannot pre-reserve one through its GitHub integration. This version is the same analysis with
the version and concept DOIs printed in Data availability and Code availability.

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
