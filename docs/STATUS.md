# Status

Written 2026-09-03. This file lives on `working-notes` because it is a development record, not
part of the release. `main` and `publication` carry the release.

## Where things stand

**The study is finished and the paper is drafted. Stop analysing.** Three council rounds and an
adversarial statistician found no surviving fatal flaw. Every remaining item is administrative
or belongs to the author.

- Repo public at `https://github.com/NirmalKumar31/rbp-protocol-calibration`.
  A **clean clone passes 649/649** verifier assertions. 648 tests. CI green on `main` and
  `publication`.
- `manuscript/paper.pdf`, 22 pages, 212-word abstract, 6 figures, 6 tables, 20 verified
  references. Build with `cd manuscript && ./build.sh`.
- Spend ~$35 of the $40 cap.

## The one outstanding item

**Zenodo DOI.** Procedure in `ZENODO.md`. Two things people get wrong: Zenodo only archives
releases created *after* the repository is switched on, and the **concept** DOI rather than the
version DOI belongs in the paper. Then edit Code availability in `manuscript/paper.tex` and
`manuscript/sections/data-availability.tex`, and rebuild.

After that the preprint goes to **bioRxiv, not arXiv**: arXiv requires an endorsement for a
first-time q-bio submitter, which is not obtainable here, and bioRxiv requires none.

## Branch layout

| branch | contents |
|---|---|
| `main` | the release; development history separated out |
| `publication` | identical to `main`; leave it alone once submitted |
| `working-notes` | this file, the old `docs/52-65`, `AGENT-CONTEXT.md`, the claims ledger `60-the-paper.md`, the review records, and the eight `manuscript/0*.md` drafts |

Nothing was deleted. If you want the reasoning behind a claim, a withdrawal, or an editorial
decision, it is on this branch.

## Two things to know before touching the code

**Run `scripts/ci_local.sh` before pushing.** CI failed four times on pushes that passed
locally, every time because the two ran different commands. That script runs exactly the
workflow's steps plus the verifier.

**The harness checks values, not provenance.** Both defects found after the assertion count
passed 600 were of that kind: 20 dinucleotide-arm datasets whose retained neural scores came
from a partition that was not chromosome-grouped, and an asymmetric convergence guard in
`lr_test` that reported a diverged fit as overwhelming evidence. Neither was caught by any
assertion. `lr_p` in `rehearsal_binding_*.csv` is unstable and unused; do not trust it.

## What is deliberately not being done

- The transform sweep, the baseline attribution and the recommendation test are reported for
  all three model classes; the order-3 analysis is 4-mer only on 30 datasets. Stated in
  Limitations.
- Extending to a third cell line or another assay. Out of scope and out of budget.
- Cutting the paper further. An editor-role review recommended six main sections; the current
  six Results subsections already match that.
