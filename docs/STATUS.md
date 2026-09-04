# Where this stands

Last updated 2026-09-04, after three review rounds.

## The paper is ready to post

28 pages, 7 tables, 6 figures, 26 references, 255-word abstract. A clean clone passes
**696/696** verifier assertions; CI green; 0 orphan numbers over 272 checked values.
`main` and `publication` are identical. Build with `cd manuscript && ./build.sh`.

## The only thing waiting on a human

**The Zenodo DOI.** `docs/ZENODO.md` has the procedure. A commented two-line sentence at the
end of `manuscript/sections/data-availability.tex` is the placeholder: uncomment it, insert the
**concept** DOI (not the version DOI), rebuild. That is the only manuscript edit left.

Later repo changes are fine: each GitHub release mints a new version under the same concept
DOI, so the paper never needs re-editing. You cannot swap files inside a published version, and
you cannot self-delete a record, so get authorship and licence metadata right the first time.

Then bioRxiv, not arXiv: arXiv needs an endorsement a first-time q-bio submitter cannot get.

## One agreed-necessary item is unpaid

**Retrain the 20 leaky dinucleotide-arm datasets.** For those 20 the committed CNN and
SpliceBERT scores came from a partition that is not chromosome-grouped. Two independent external
reviews and an internal adversarial referee all demand a retrain rather than a disclosure.

- 200 fold-runs (20 datasets x 5 folds x 2 models)
- those datasets hold **37% of the dinucleotide arm's rows**, so ~**$7** on Modal, not the ~$4
  a run-count estimate suggests
- the defect is **conservative**: dropping the 20 WIDENS every span (5.42 to 6.02, 7.63 to 8.33,
  3.76 to 3.92) and every restricted value sits inside the full-panel interval

So the paper is postable without it. It just leaves the most quotable reviewer objection
standing. The full panel is primary and the 74 are the sensitivity, which is stated in both
Results and Limitations.

## Free work that was offered and not run

All local CPU, no cloud:

| item | effort | why |
|---|---|---|
| gene/transcript-clustered CV | ~3 h | the only one that could ADD a result, not just close a hole |
| train-fold-only standardisation | ~1 h | removes a certain reviewer comment; currently label-free but improper |
| positive-set intersection as a sensitivity | ~1 h | Jaccard is 0.9972, so this is a sensitivity, not a primary panel |
| matching-algorithm robustness | ~2 h | fair ask given the paper's own thesis |
| study-design figure | ~1 h | the methods take 8 pages to grasp and should take one figure |

## Two things not to get wrong

**Read the rendered PDF, not the source.** Three figure captions described panels their figures
did not contain, and no source-level check can see that. Extract with `pypdf` and read it.

**Six claims were wrong and are fixed.** The worst: the composition block has rank **18**, not
15, because dinucleotide counts sum to L-1 and recover base counts only up to the terminal base.
That sentence was the stated reason the primary contrast's direction is implied by the design.
The conclusion survives on the asymmetry argument instead. The others are the external benchmark
being external in construction only, the withdrawn capacity claim, headroom normalisation not
being a monotone transformation, undefined "apparent AUROC", and a backwards variance argument.

## Branches

- `main`, `publication`: identical. Leave `publication` alone once submitted.
- `working-notes`: this file, the development record, council records, draft markdown.
