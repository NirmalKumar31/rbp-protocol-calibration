# Today's in-depth final audit

**Repository:** `rbp-repro` / `NirmalKumar31/rbp-protocol-calibration`  
**Audited commit:** `3b5a2d3e6f944d2d397ecee970983b9c9b022345`  
**Audit date:** 2026-09-06 (America/Denver)  
**Mode:** Read-only review. No research code, manuscript text, data, tables, figures, configuration, or release metadata was changed. This report is the only file added by the auditor.

## Final verdict

**No: the exact repository state at the audited commit is not ready to publish or tag.**

The scientific study itself is publishable. Its central, properly narrowed conclusion is well supported: in this ENCODE eCLIP panel, the *measured* incremental contribution of a sequence model is strongly indexed to negative-set construction. The paper is unusually candid about estimator bias, the failed external directional replication, exploratory neural results, and the limits of generalisation. I found no evidence that the headline 4-mer result is numerically false or that the study must run the unaffordable neural cross-fitting before a preprint.

However, Claude's statement that the audit list is closed except for two decisions is factually wrong. The latest pushed commit has a **red CI run**, its new cache-idempotence release gate fails in a clean archive, and several manuscript/release defects remain. Publishing this exact commit would knowingly publish a release whose own newly advertised reproducibility check fails.

### Readiness scores

These percentages are judgement calls, not statistical probabilities:

| Dimension | Current readiness | Assessment |
|---|---:|---|
| Scientific question and value | **92%** | Clear, useful calibration question with practical implications for RBP benchmarking. |
| Evidence for the narrow primary claim | **90%** | Strong paired design, primary cross-fitted 4-mer result, clustered uncertainty, multiple robustness checks, and an external protocol-dependence replication. |
| Ability to survive scientific scrutiny | **85–90%** | Likely to survive if claims stay narrow; reviewers can fairly challenge conditional uncertainty, model scope, panel representativeness, and incomplete neural cross-fitting. |
| Reproducibility engineering | **82%** | Excellent released evidence and verification coverage, but the newest reproducibility gate currently fails and has blind spots. |
| Manuscript accuracy and internal consistency | **80%** | The main argument is coherent, but several direct wording contradictions remain. |
| bioRxiv package, at this exact commit | **70%** | Scientifically postable after the blockers below; not safe to upload now. |
| Journal submission package | **60–65%** | A 58-page, 440-word-abstract manuscript still needs venue-specific restructuring and substantial editing. |

After the P0 and P1 items below are closed, I would regard it as ready for **bioRxiv**. That is different from saying it is already formatted and edited for a particular journal.

## P0 — must fix before any tag, DOI deposit, release, or preprint upload

### P0.1 The current commit fails its own cache-idempotence gate

On a fresh `git archive` of `3b5a2d3`, I ran:

```bash
python scripts/cache_idempotence.py
```

All 33 entry points executed, but the command exited 1:

```text
33 of 33 entry points ran, tolerance 1e-09

2 ENTRY POINT(S) DO NOT REPRODUCE THEIR COMMITTED TABLE:

gene_clustered_cv
  gene_clustered_cv.csv: NOTE changed on three carried rows
window_centring
  window_centring.csv: NOTE changed on two carried rows
```

GitHub Actions independently failed at exactly the same step in [run 34086455394](https://github.com/NirmalKumar31/rbp-protocol-calibration/actions/runs/34086455394). The full test suite, lint job, manuscript job, and prior steps passed, but the overall run is red.

The immediate cause is `src/rbp/utils/carry.py:52–54`: every cache run appends `"; not recomputed here; the input is absent"` to the already-carried note. A second run therefore necessarily differs from the committed table. Blank notes also arrive through pandas as `NaN`; `str(row.get("note") or "")` does not reliably normalise that case because `NaN` is truthy.

This directly contradicts the latest commit message, which asserts that all 33 entry points reproduce cleanly. Do not tag this commit.

**Acceptance criterion:** run the gate twice in succession in a new `git archive`; both runs must exit 0, and the archive's tables must be restored exactly after each run.

### P0.2 The new gate does not yet establish everything its documentation claims

Even after fixing the repeated suffix, `scripts/cache_idempotence.py` has material false-pass paths:

- Lines 143–146 classify a non-zero entry-point exit as **skipped**, not failed. A broken script can therefore be omitted from the evidence and the overall gate can still pass.
- Lines 86–117 compare the set of `check` rows, `note`, and `value`, but not the complete schema or other claim-bearing fields such as `n`, `ci_low`, `ci_high`, or analysis-specific columns. The docstring calls the comparison structural and numeric more generally than the implementation supports.
- Lines 148–159 examine only top-level `*.csv`; they ignore TSV files and nested table paths.
- Restoration copies snapshot files back but does not remove new files created by an entry point. Such a file can persist, contaminate the next entry point, and remain after the gate. The comment that the script “restores what it touched” is therefore too strong.
- There are no focused unit tests for `rbp.utils.carry` or the cache-idempotence comparator. The commit message says the gate was attacked manually, but that regression is not preserved in the test suite.

**Required fix:** non-zero exits must fail; compare exact column sets, row identity/uniqueness, text fields, and every numeric field with explicit tolerance; cover all advertised output paths; restore by replacing the whole snapshot scope, including deleting additions; add tests that deliberately exercise each failure path.

### P0.3 The local CI mirror omits the new failing gate and masks another drift

`.github/workflows/ci.yml` says changes must be mirrored in `scripts/ci_local.sh`. The workflow added `python scripts/cache_idempotence.py`, but `scripts/ci_local.sh` does not run it. Thus the command advertised as locally reproducing CI would pass while CI fails.

In addition, local CI runs `scripts/column_dictionary.py` **without** `--check`. That rewrites the dictionary instead of detecting a stale committed dictionary, whereas GitHub CI correctly uses `--check`.

**Acceptance criterion:** the local script must invoke the same new gate and use `column_dictionary.py --check`. A clean local run and the GitHub run for the final commit must both be green.

### P0.4 Resolve the AI-use disclosure before submission

`manuscript/paper.tex:201–222` explicitly labels this an open item. The rendered disclosure names only Anthropic Claude for coding and drafting, while the source comment records six external generative-system audits that materially changed code, tables, and manuscript statements.

This is not merely a stylistic choice. The author must inventory actual use, check the target venue's current policy, and make the disclosure complete and accurate. The sentence “no reference was generated by a language model” is a personal factual warranty and should remain only after the author checks every bibliography entry.

No automated assistant can make this decision on the author's behalf. The release must not retain an “OPEN ITEM ... BEFORE SUBMISSION” marker while being described as final.

### P0.5 Resolve version, tag, and archival DOI consistently

There is no Git tag or project DOI. `pyproject.toml` and `CITATION.cff` both say `0.9.0`. A publication snapshot needs an immutable identifier, and the version in both files must match the Git tag.

The DOI instructions contradict one another:

- `SUBMISSION.md:157–161` says to put the **concept DOI only** in the manuscript and calls it the only necessary edit.
- `docs/ZENODO.md:81–85` correctly says the **version DOI** identifies the exact analysed snapshot and recommends giving both.
- `docs/ZENODO.md:93–94` then switches back to recommending the concept DOI in the manuscript.
- `manuscript/sections/data-availability.tex:38–45` says to insert **both**.

For research reproducibility, cite the version DOI and matching tag for the exact snapshot, with the concept DOI additionally for discovery. Choose and document one non-contradictory release sequence. If using automatic GitHub–Zenodo minting, acknowledge the DOI circularity and create the follow-up DOI-bearing release as documented; alternatively reserve a DOI through a manual Zenodo draft before the final archive.

I agree that `v1.0.0` is a sensible version for the first archival publication, but only after the blockers are fixed. Do not use a tag to make a failing commit look final.

### P0.6 Make an explicit decision about the historical billing-account identifier

The history scan passes because the known finding is stable, not because the history is clean. It reports four commits, nine commit/file pairs, and eleven diff lines containing a live GCP billing-account ID. `SECURITY.md` accurately explains that this is an identifier rather than a credential, but the owner must choose one of the documented options before making an immutable DOI snapshot.

This is an owner decision, but it remains a release blocker until recorded. If retaining history, explicitly accept the exposure and confirm the account is safe to name publicly. The cleanest security outcome is to move future work to a different billing account before archival release.

### P0.7 The release archive omits the licence file it promises for `results/`

`results/LICENSE` exists only as an ignored local file. It is not tracked and is absent from `git archive HEAD` because `.gitignore` excludes `/results/*` without re-including the licence. The root `LICENSE` explicitly promises that both `results/` and `data/evidence/` carry their own licence so either directory remains licensed when extracted alone. Only `data/evidence/LICENSE` is committed.

The root licence still declares results to be CC BY 4.0, so the legal intent is visible at repository level. The defect is that the release artefact contradicts its explicit packaging claim.

**Required fix:** re-include and commit `results/LICENSE`, then verify it appears in `git archive` and in the Zenodo payload.

## P1 — fix before calling the package publication-ready

### P1.1 Correct stale “relation in the title” language

The title was deliberately changed to the robust claim, “Apparent sequence-model contribution depends strongly on negative-set construction.” The Discussion still says at `manuscript/sections/discussion.tex:159–160`:

> The relation in the title ... does not hold on the one externally constructed benchmark.

That is now false. Protocol dependence *does* replicate on the 135-dataset external complement (1.69-fold, 95% CI 1.41–2.03), while the **directional inverse relationship between apparent AUROC and contribution** does not. The current sentence conflicts with the preceding paragraph and the revised title.

Replace “relation in the title” with the exact directional relation that failed. This is a substantive clarity correction, not a cosmetic edit.

### P1.2 Resolve the external model-only-AUROC contradiction

`manuscript/sections/results.tex:1356–1363` says Horlacher's released sets carry no model-only AUROC, that difficulty “can only” be read from the composition baseline, and that the 4-mer's own AUROC is unavailable. Lines 1393–1403 then report that this study measured the 4-mer's own AUROC on those same sets (0.7647 and 0.8594).

The intended distinction is understandable: the external release did not supply model-only AUROC, whereas this study subsequently computed it. The current wording nevertheless reads as a direct contradiction. Say explicitly that it was **not supplied by the release** and that the study's separately computed values are presented below. Update the table caption as well.

### P1.3 Correct the baseline feature count in the supplement

`manuscript/supplementary.tex:51–53` says the baseline has 19 features and then defines it as four mononucleotide frequencies, sixteen dinucleotide frequencies, and entropy—21 columns.

The implementation and Methods are clear: one level is dropped from each frequency family, leaving **3 + 15 + 1 = 19** columns. Correct the supplement to match `manuscript/sections/methods.tex:221–236` and `src/rbp/eval/nested.py:117–141`.

### P1.4 Fix the ambiguous baseline-indexing statement

`manuscript/sections/discussion.tex:150–151` says:

> Every magnitude reported here is indexed to a mono- and dinucleotide baseline. The protocol dependence is not so indexed.

The reported nested contributions and their spans are, by definition, indexed to that baseline; the baseline-order sensitivity demonstrates this. If the intended claim is that the **qualitative existence** of protocol dependence is robust to alternative baselines/estimands, say precisely that. As written, the second sentence contradicts the first and overstates invariance.

### P1.5 Make the main and supplement titles identical

The main manuscript title begins “Apparent sequence-model contribution depends strongly...”, but `manuscript/supplementary.tex:33–35` uses the older title “What a sequence model appears to contribute is set by...”. The supplement metadata uses a third shortened form. A submission system and readers should be able to associate the files unambiguously. Use the exact final article title in both visible titles and metadata.

### P1.6 Repair the two supplement overfull boxes and normalize page size

A clean build succeeds, but `supplementary.log` contains:

- `Overfull \hbox (20.81735pt too wide)` at lines 73–90: the Figure/Table mapping table.
- `Overfull \hbox (21.8744pt too wide)` at lines 339–343: the long supplementary CSV/path text.

The main manuscript has one non-fatal underfull box at `results.tex:452–466`. Visual review of all pages found no clipped figures, overlaps, missing glyphs, or accidental blank pages, but the supplement's overflows are visible and should be corrected before upload. `build.sh` checks undefined citations/references but does not fail on overfull boxes, so its successful message is not a layout-clean certificate.

The main PDF is US Letter (612 × 792 pt); the supplement is A4 (595.28 × 841.89 pt). Pick one page size, ideally the target venue's requirement, and use it consistently. Add supplement PDF keywords while normalising metadata.

### P1.7 Make provenance wording match behaviour

`scripts/provenance.py --help` says `--check` fails if “any table is unattributed,” but the command exits 0 with two unattributed files:

- `results/tables/unattributed/backend_replication.csv`
- `results/tables/unattributed/substitution_baseline.csv`

Their README says they were moved out of `results/tables/`, although they remain under its `unattributed/` subdirectory and are included in the manifest. They are not used by this paper and their retention is transparently explained. The problem is the mismatch among help text, directory prose, and actual acceptance behaviour.

Either make unattributed files fail as advertised, or explicitly document this quarantined exception everywhere. For a paper-specific archival release, consider moving the earlier-study artefacts to a separately versioned archive; do not delete evidence silently.

### P1.8 Correct stale release documentation

The `.gitignore` comment describes `data/evidence/` as 1,120 files and 35 MB and says 95/95 deep-model datasets reproduce. The current tree contains roughly 3,496 tracked evidence files and 71 MB; the paper's three-arm panel is 94 datasets, with a 95th candidate excluded from full three-arm results. Comments do not affect execution, but this repository repeatedly presents comments as audit history, so stale facts erode confidence.

Review `SUBMISSION.md` after all fixes. Its statement that DOI insertion is the “only manuscript edit required” is already disproved by the manuscript issues above.

## P2 — scientific and editorial limitations to preserve, not hide

These are not reasons to block a carefully framed preprint. They are the points a serious reviewer is likely to probe, and the manuscript should continue to state them plainly.

### Scientific scope

- The primary cross-fitted estimator is available for the k-mer models, not the CNN or SpliceBERT. Neural spans use the biased two-stage estimator and must remain explicitly exploratory.
- Neural initialisation was not fully seeded across the original sweeps, and draw/initialisation variability is not propagated through every reported neural interval.
- The protein-clustered bootstrap correctly handles proteins observed in both cell lines, but intervals remain conditional on the selected panel, one fold partition, chosen hyperparameters, retained positives, and mostly one negative draw.
- The panel is a deterministic systematic subset by dataset size, not a probability sample. Generalisation should remain limited to the analysed ENCODE eCLIP K562/HepG2 setting.
- Exact 32-mer overlap deletion is a useful leakage sensitivity, not a tunable sequence-homology clustering. Diverged homology can remain.
- The external held-out complement is valuable, but the benchmark family and overlapping subset were already known before the disjoint complement was scored. Calling it externally constructed is fair; calling it a fully prospective, independent validation would not be.
- The external test replicates protocol dependence but not the paper's internal inverse direction. This failure is a strength of the reporting and must remain attached to any abstract, press, or README claim.
- The seven-source negative-set survey is targeted and hand-selected, not systematic. “None of seven surveyed” is supported by the committed classifications; “the field never reports this” would be too broad.
- The class-ratio sensitivity and model-capacity ladder were specified but not run. They are not needed for the narrow primary claim, provided the manuscript does not revive a capacity-ordering claim or imply class-ratio robustness.
- A full raw rerun does not follow one uniform command for all arms: `run.sh all` covers the dinucleotide arm, while GC and bias-aware neural sweeps have separate cloud/Modal drivers. Released per-window scores allow exact downstream recomputation, but this is evidence-based reproducibility rather than a complete low-cost raw replication.

### Statistical interpretation

- The central quantities are descriptive protocol contrasts, not causal effects. The manuscript generally gets this right.
- Ratios of small positive increments can look dramatic. The paper responsibly reports additive contrasts, alternative estimands, and confidence intervals; retain those adjacent to fold-change claims.
- DeLong testing is invalid for the nested-model null addressed here. The manuscript appropriately treats dataset-level p-value counts as descriptive and bases headline inference on clustered bootstrap intervals; do not drift back to “significant in X datasets” as primary proof.
- Numerous robustness and exploratory analyses are not multiplicity-adjusted. Keep them labelled exploratory and avoid treating the count of successful sensitivities as a formal family-wise test.
- The manuscript-number audit is a regression/traceability tool, not proof of correctness. It reports high collision/false-negative rates: 28.9% occupancy of the four-decimal grid, 96.0% of the three-decimal grid, and 37.0% of integer values from 10–1000. Its zero-orphan result is useful, but should never be described as validating the scientific meaning of every number.

### Writing and journal suitability

The prose is technically precise but too dense for a final journal submission. Automated prose review found approximately 22,834 section words, 100 sentences of at least 45 words, several very long paragraphs, and a Results section containing more than half the body text. The 58-page paper and 440-word abstract may be acceptable to bioRxiv, but “bioRxiv has no limit” does not mean the manuscript is optimally readable or journal-ready.

For a journal version, reduce repetition among Results, Discussion, comments, and captions; shorten sentences; move secondary robustness detail to the supplement; and target the selected journal's article structure and limits. This is venue-dependent, but the need for substantive editing is not merely optional polish.

## What passed

The negative verdict should not obscure substantial strengths. At the audited commit:

- Full GitHub unit suite passed: **796 collected tests** in the full-suite job.
- Local CPU-compatible unit suite passed with six expected skips and the two torch-dependent modules excluded.
- `ruff check .` passed.
- Every tracked shell script passed `bash -n` in CI.
- Terraform configuration validated successfully.
- Repository packaging succeeded and import-tested nine modules.
- `scripts/verify.py --local results/tables` passed **1083/1083** assertions: 945 for this paper, 136 for the earlier variant study, and two harness self-checks.
- `scripts/audit_manuscript.py` checked 690 decimal values and 500 integers and reported zero orphans, with its limitations disclosed.
- The raw-input manifest checks 251 objects with size and MD5 and rejects unsafe paths/secrets.
- The provenance manifest matches 158 committed tables: 25 raw-reproducible, 31 evidence-recomputable, 33 frozen-cache, 44 frozen-only, 23 cloud-produced, and two quarantined unattributed files.
- The column dictionary matches all 158 tables and 1,607 columns.
- CFF 1.2 metadata validates, and `pyproject.toml`/`CITATION.cff` currently agree at 0.9.0.
- The paper builds to 58 pages and the supplement to 12 pages, with no undefined citations or references.
- Visual inspection found the seven main figures and ten supplementary figures readable, with no clipping or broken rendering.
- The newest 2026 Tourne et al. citation and DOI metadata were independently checked against the DOI record and agree on title, authors, journal, volume, issue, article number, and year.
- The history scanner found no private keys, service-account key material, AWS keys, GitHub tokens, Slack tokens, Google API keys, or OAuth tokens. The billing identifier is the disclosed exception.

## Scientific bottom line

The study makes sense and is worth publishing. Its likely impact is methodological rather than biological: it demonstrates that a commonly compared “model contribution” is partly a property of benchmark construction, quantifies the effect, shows an estimator floor, and gives an actionable reporting/cross-fitting recommendation. This is a useful result for computational biology benchmarking.

It will survive reasonable scrutiny **if and only if** the claim remains narrow:

> In this analysed eCLIP setting, measured nested contribution is strongly protocol-indexed; the primary cross-fitted 4-mer span is robust, while neural magnitudes are exploratory and the internal inverse directional relationship fails on the external benchmark.

It should not claim universal causal effects of negative sampling, a universal inverse relation, a model-capacity law, fully cross-fitted evidence across all model classes, or generalisation beyond the analysed assays and panels.

No finite audit can make science “100%” certain. Publication readiness means there are no known release-breaking contradictions, all remaining limitations are accurately disclosed, and independent gates pass. This commit does not yet meet that standard, but it is close and the remaining technical/scientific-writing fixes do not require a new expensive experiment.

## Exact final checklist

Do not publish until every item below is checked:

- [ ] Fix `carry.py` so carried notes are idempotent, including pandas `NaN` handling.
- [ ] Harden `cache_idempotence.py` against non-zero exits, incomplete-field comparisons, nested/TSV outputs, duplicate keys, and incomplete restoration.
- [ ] Add focused tests for carry/idempotence failure modes.
- [ ] Add the gate to `scripts/ci_local.sh`; change column dictionary invocation to `--check`.
- [ ] Run the cache gate twice in a clean archive; both passes must leave no table differences.
- [ ] Obtain a fully green GitHub Actions run on the exact final SHA.
- [ ] Correct “relation in the title,” external model-AUROC availability, baseline feature count, and baseline-indexing wording.
- [ ] Align supplement title/metadata with the paper.
- [ ] Remove the two supplement overfull boxes and normalise paper size.
- [ ] Reconcile provenance help/README behaviour for quarantined unattributed files.
- [ ] Commit `results/LICENSE` and verify it appears in the release archive.
- [ ] Refresh stale `.gitignore` evidence counts.
- [ ] Author approves a complete, venue-compliant AI-use disclosure and verifies the bibliography warranty.
- [ ] Author records the billing-ID decision.
- [ ] Choose version; make `pyproject.toml`, `CITATION.cff`, and tag agree.
- [ ] Reconcile DOI instructions; mint and cite the exact version DOI plus concept DOI as decided.
- [ ] Build both PDFs from the final clean archive and inspect logs and rendered pages.
- [ ] Confirm the final archive contains PDFs, supplementary CSV, per-directory licences, citation metadata, tables, figures, and released evidence.

## Message to send Claude

Copy and send the following verbatim:

> Do not stop and do not tag or create a DOI yet. Audit commit `3b5a2d3e6f944d2d397ecee970983b9c9b022345` is not publication-ready. Its GitHub Actions run 34086455394 is red, and `python scripts/cache_idempotence.py` fails in a clean `git archive` on `gene_clustered_cv` and `window_centring` because `rbp.utils.carry.carried()` appends the missing-input suffix again on every run. Fix this class, not only the five rows.
>
> Harden the gate before trusting it: a non-zero child exit must fail, not be reported as skipped; compare exact schemas, unique row keys, all text fields, and all numeric fields (`value`, `n`, confidence bounds, and analysis-specific numeric columns) with documented tolerance; cover every advertised CSV/TSV/nested output; and restore the complete output tree including deleting new files. Add focused regression tests for repeated carry, blank/NaN notes, child failure, dropped/added/duplicate rows, changed `n`/CI, new files, and restoration after failure.
>
> Update `scripts/ci_local.sh` to run the new gate and use `column_dictionary.py --check`, so it actually mirrors CI. Prove the result by running the cache gate twice in succession in a fresh archive and showing that both runs exit 0 and leave no differences. Then obtain a green GitHub Actions run for the exact final SHA.
>
> Fix the remaining manuscript contradictions: (1) Discussion's “relation in the title ... does not hold” is stale because the revised title's protocol-dependence claim *does* replicate externally; name the inverse directional relation instead. (2) Results first says model-only AUROC is unavailable/can only be inferred from composition and later reports the study's computed 4-mer AUROCs; distinguish “not supplied by Horlacher” from “computed here,” including the caption. (3) Supplement says 19 features but lists 4 + 16 + entropy; state the actual 3 + 15 + entropy after dropping one level per family. (4) Clarify “protocol dependence is not indexed” because all reported nested magnitudes are baseline-indexed. (5) Make the supplement's visible title and metadata exactly match the paper.
>
> Rebuild and remove the supplement's two overfull boxes (20.81735 pt at lines 73–90 and 21.8744 pt at 339–343), make paper and supplement use one page size, and add consistent metadata. Do a rendered-page inspection, not only a successful LaTeX exit.
>
> Fix release packaging: `results/LICENSE` exists locally but is ignored and absent from Git/`git archive`, despite the root licence promising it; re-include and commit it. Reconcile provenance documentation: `--help` says unattributed tables fail, but two quarantined files pass, and their README incorrectly says they are outside `results/tables`. Refresh the stale `data/evidence` count/size comment in `.gitignore`.
>
> Reconcile all DOI guidance. The exact analysed snapshot needs its version DOI and matching tag; the concept DOI can be given additionally for discovery. `SUBMISSION.md`, `docs/ZENODO.md`, the manuscript comment, and `CITATION.cff` must tell one consistent workflow. Recommend `v1.0.0`, but only after the author settles the AI-disclosure and historical billing-ID decisions. Ensure `pyproject.toml`, `CITATION.cff`, and the tag agree.
>
> Do not run the ~$573 neural cross-fitting or invent new claims. The narrow primary science is publishable without it; keep CNN/SpliceBERT explicitly exploratory and preserve the failed external directional replication and all conditional-inference limitations. Do not describe class-ratio or capacity-ladder robustness as tested.
>
> Return one final evidence packet: final commit SHA; clean `git status`; exact output of two cache-idempotence runs from a clean archive; full test count; 1083/1083 verifier result; provenance/column/raw-input checks; CFF validation; Terraform validation; PDF page counts and zero undefined/overfull warnings; archive inventory proving both licences and supplements are present; and the green GitHub Actions URL. Separately list the exact author decisions still required. Do not say “everything is done” unless every mechanical item above has evidence and the only remaining items genuinely require the author's decision.

## Release decision

**Decision: HOLD.** Fix the P0/P1 defects, obtain a green final SHA, settle the owner decisions, and then release. No new paid scientific run is required for the paper's narrow primary claim.
