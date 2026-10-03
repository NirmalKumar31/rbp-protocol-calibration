# Changelog

## 1.0.9, the manuscript matches the paper of record

No code changed and no published number moved. The numeric token count of the rebuilt PDF is
unchanged at 2652.

A claim-by-claim comparison of the posted Research Square paper against this repository found
**no numeric drift**: 24 published values, from the abstract and Tables 1, 2, 15 and 16 and the
external span, all match the committed tables exactly. One textual inconsistency was real, and
this release carries its correction into the archive.

- **The declarations now match the posted paper.** It runs Funding, Competing interests,
  Ethics. The repository ran Funding, Ethics, Conflicts of Interest, and declared "no conflicts
  of interest" where the posted paper declares "no competing interests".
- **The cause is recorded beside the heading.** That ordering came from 1.0.3, whose entry
  below says it converged the repository on the posted paper. It converged on Preprints.org.
  Once Research Square became the version of record, the repository was matching an abandoned
  submission target against the paper of record.
- **Documentation gained estimator labels.** Five places stated a span without naming which
  estimator produced it: the README summary table, `docs/PANELS.md`, `PORTFOLIO.md` and two in
  the dashboard. Every value was correct and traced to a table; what was missing was whether
  it came from the two-stage or the cross-fitted estimator. Three readers divided a two-stage
  figure, compared it to the cross-fitted 4.84 headline and reported a contradiction. No gate
  could catch this, because both numbers are real, which is the limitation the paper's own
  Limitations section states.

**Data availability still cites v1.0.0** as the analysed snapshot. That is correct rather than
stale: v1.0.0 is the state that produced the numbers, and the concept DOI resolves to the
latest version, which is why the paper cites the concept DOI and not a version DOI.

## 1.0.8, the audit stops reading URLs as claims

No code that produces a result changed, and no published number moved.

A video embed added to `README.md` broke CI. The asset link
`github.com/user-attachments/assets/d38e5be5-...-8846-...` carries `8846` inside its UUID;
`scripts/audit_manuscript.py` extracted that as an integer claim, could trace it to no table,
and failed the gate. The README was correct and the audit was wrong.

- **A number inside a URL is not a claim.** No reader reads a path segment as a quantity.
  `IDENTIFIER` already skipped DOIs, Zenodo ids, ENCODE accessions and version strings by
  inspecting surrounding characters, but it cannot enumerate every opaque identifier that will
  ever appear in a link. Whole URLs are now blanked before scanning, which is the general form
  of the same rule, substituted with spaces so reported column positions stay honest.
- **The gate is not weakened**, and that was checked rather than assumed. Integers checked fall
  from 494 to 493, which is the one URL digit and nothing else. A fabricated value and a
  fabricated count injected into README prose are both still caught; a URL full of digits is
  ignored.

This release also carries the README video link, which was added after the 1.0.7 tag and so
was not in that archive.

## 1.0.7, one preprint named

Metadata only. No code, no results, no manuscript text changes.

Research Square, `10.21203/rs.3.rs-10988414/v1`, is the only preprint this repository names.
The `CITATION.cff` note describing an earlier posting is removed, and the 1.0.3 and 1.0.4
entries below are now venue-neutral where they previously named a server and a DOI.

**This rewrites changelog history, which is normally not done here.** It is recorded rather
than done quietly: those entries described what happened on the day they were written. The
earlier record is not affected by any of this. A posted DOI cannot be unposted, and the
archived snapshots for 1.0.3 through 1.0.6 are immutable and still contain the original text.
What changes is what this repository points a reader towards from here.

## 1.0.6, Research Square is the version of record

Metadata only. No code, no results, no manuscript text changes.

The paper is posted on Research Square: DOI `10.21203/rs.3.rs-10988414/v1`, 14 September 2026,
CC BY 4.0, published by Springer Nature. That is now the version this repository cites.

- **`CITATION.cff`** `preferred-citation` carries the Research Square DOI and URL.
- **`.zenodo.json`** declares it as `isSupplementTo`, so the archive points at the version of
  record rather than at the earlier one.
- **`README.md`** links it as the paper, and as the DOI to cite for the study.

## 1.0.5, the lint fix 1.0.4 needed

1.0.4 was tagged and released with a failing `lint` job. The drift-guard test added in that
release assigned a lambda to a name, which ruff rejects as E731. The other four CI jobs
(`test`, `full-suite`, `manuscript`, `docker-cpu`) passed, and no result, number, figure or
table was affected: the rule is a style rule and the file is a test.

- **`collapse` is a `def` rather than a lambda** in
  `test_zenodo_json_and_citation_cff_cannot_drift`. `ruff check .` passes repository-wide.
- **The 1.0.4 tag and its Zenodo archive are left exactly as they are.** Moving a tag that a
  DOI has already been minted against would break the one correspondence a version DOI
  exists to guarantee, which is a worse defect than the one being fixed. `10.5281/zenodo.22712185`
  therefore archives a tree whose lint job fails, and this entry is the record of why.
- **Cause, recorded because it will recur otherwise:** 1.0.4 was verified with `pytest` and
  `scripts/release_consistency.py` run separately, neither of which lints. `scripts/ci_local.sh`
  reproduces all five CI jobs including ruff and would have caught it before the tag existed.
  It is the only pre-tag check worth running.

Verified before tagging this time: all CI-local checks pass, 37 of 37 cached entry points
reproduce their committed tables at 1e-09, and both tracked PDFs are byte-identical to a
clean build.

## 1.0.4, the archive points at the paper

Metadata only. No code, no results, no manuscript text changes.

The Zenodo record carried one related identifier, a link to the GitHub tree, and nothing
saying which paper the deposit belongs to. `CITATION.cff` gained the preprint DOI in 1.0.3,
but Zenodo has no field that maps a CFF `preferred-citation` onto a related identifier, so
the archive stayed silent about the paper it supports.

- **`.zenodo.json` added**, declaring the posted preprint as `isSupplementTo`, resource type
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
paper say the same thing, which they did not between 1.0.2 and the preprint posting.

- **The paper is posted**, CC BY 4.0. `CITATION.cff` records it in
  `preferred-citation`, which previously carried no
  DOI at all and so left no machine-readable trace of where the work appeared.
- **Use of AI tools.** Methods gains a subsection naming the tools, what they were used for
  and who is responsible. The preprint server requires the declaration; it belongs in the
  paper regardless, and the repository copy lacked it while the submitted copy had it.
- **Conflicts of Interest.** The declaration was titled "Competing interests" and sat before
  Ethics. It is now titled and positioned as the submitted version has it.
- **A false comment removed.** `CITATION.cff` asserted that the software title and the paper
  title "carry the same title on purpose". They differ, deliberately, and have since this
  record became `type: software`. The comment described a decision that had been reversed:
  a pointer that read confidently and was wrong, which is the defect class this paper is
  about.

`manuscript/paper.tex` and `manuscript/sections/methods.tex` are now byte-identical to the
source that was submitted. Every published number is unchanged.

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
- The analysis code, the cloud definitions that ran it, and the container build definitions.
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
