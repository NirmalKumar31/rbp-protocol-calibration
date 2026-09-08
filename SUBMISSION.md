# Submission package

Everything needed to post the preprint, and everything needed to verify every number in it
offline against committed evidence. That is not the same as rebuilding every number from raw
inputs. `results/tables/PROVENANCE.csv` classifies every released table into
raw-reproducible, evidence-recomputable, frozen-cache, frozen-only, cloud-produced and
unattributed, and the counts are in that file rather than retyped here, because retyped counts
in this file are exactly what an external review found eight of. Most tables are verified
rather than reconstructed. Use "verified against released evidence" and "reconstructed from raw
inputs" as distinct claims throughout; this file said "reproduce every number" and meant the
first. Nothing here is a summary of the science, which is in `manuscript/paper.pdf`.

## What to upload

| bioRxiv asks for | file |
|---|---|
| Manuscript PDF | `manuscript/paper.pdf` (63 pages) |
| Abstract (paste into the form) | abstract of `manuscript/paper.tex`, 536 words, no markup |
| Supplementary tables | `results/tables/supplementary_table_s1.csv` and the per-dataset tables listed below |
| Source, if requested | `manuscript/` is self-contained: `paper.tex`, `sections/`, `figures/`, `build.sh` |

**Venue: bioRxiv.** **Subject area:** Bioinformatics. **Licence:** CC BY 4.0, matching the
licence on `results/` and `data/evidence/`. **Type:** New Results.

bioRxiv's requirements were checked against its submission guide rather than assumed, because
three reviews raised length and supplement questions that turn entirely on the venue. What the
guide actually states is that it **does not specify** page, abstract, figure, table or reference
limits, that it accepts a single PDF containing text and display items, and that supplementary
material may be supplied as separate files. That is an absence of stated limits, which is not
the same as a guarantee that none exist -- an earlier version of this paragraph asserted the
stronger claim and a review was right to object.

**Check the live submission form on the day of upload.** Nothing below depends on a limit being
absent; it depends only on nothing being misdescribed.

What that does NOT settle: a journal submission afterwards will impose its own limits, and the
ones that would bite are the abstract (536 words against a typical 250) and the main text (62
pages, 17 table environments). The long-form abstract and the guidance for cutting are kept
where they can be found again -- see the note above `\begin{abstract}` in `paper.tex` -- and
`results/tables/PROVENANCE.csv` already identifies which tables are secondary and would move to
a supplement first.

## The manuscript

`manuscript/paper.tex` builds `paper.pdf` with repeated runs of `pdflatex` until the
cross-references reach a fixpoint, which is up to four and is not reliably two, or in one
step:

```
cd manuscript && ./build.sh
```

`build.sh` derives the figure list from `sections/*.tex` and copies those figures out of
`results/figures/`, so the directory is a self-contained upload that cannot go stale against a
figure added to the text. It **fails** rather than warns on an undefined reference or citation.

Structure follows the convention of the target literature (Horlacher *et al.* 2023 and Chen
*et al.* 2024, both *Briefings in Bioinformatics*): Abstract, Introduction, Materials and
Methods, Results with declarative subsection headings, Discussion with Limitations, Conclusion,
then the declarations and 31 references.

`manuscript/` contains only what the submission needs: `paper.tex`, `sections/`, `figures/`,
`build.sh`, the built `paper.pdf` and `supplementary_table_s1.csv`, which ships with the
manuscript because a submission whose supplementary file is a repository path is not a
submission. Drafting notes and the record of editorial decisions are
on the `working-notes` branch.

## Main display items

| item | file | section |
|---|---|---|
| Figure 1 | `results/figures/f10_three_protocols.pdf` | three protocols, and the reversal |
| Figure 2 | `results/figures/f11_scale_sweep.pdf` | eight reparameterisations |
| Figure 3 | `results/figures/f12_protocol_or_baseline.pdf` | baseline attribution |
| Figure 4 | `results/figures/f9_deep_contrast.pdf` | model classes |
| Figure 5 | `results/figures/f16_order_profile.pdf` | baseline order, and where the baseline stops being one |
| Figure 6 | `results/figures/f14_external_validation.pdf` | externally constructed benchmark |
| Figure 7 | `results/figures/f15_recommendation.pdf` | the recommendation, in and out of sample |

The numbering is the order the figures appear in the text, which `manuscript/build.sh` derives
from `sections/*.tex`. Adding a figure to the text renumbers this table; it is not maintained
by hand and should be reread rather than trusted.

The manuscript's tables are typeset from the committed tables below; a count is deliberately
not given here, because it went stale twice. Every figure PDF is pure vector with embedded
TrueType fonts, which `tests/unit/test_figure_output.py` asserts rather than this file
claiming it.

## Supplementary material

Column definitions for every released table are in `results/tables/COLUMNS.csv`, one row per
file per column, generated by `scripts/column_dictionary.py`: dtype over the whole column, unit,
key flag, row and missing counts, numeric range, a definition where the release's controlled
vocabulary covers the name, and the producing script for every column either way.
`results/tables/COLUMNS_SUMMARY.csv` holds the counts, which are generated and not typed.
`results/tables/SCHEMA.md` explains the two common shapes, the units, the missing-value
conventions and the protocol-suffix mapping.

Coverage is not implied, it is counted: `COLUMNS_SUMMARY.csv` records how many columns carry
a prose definition and how many name a producing script, and the shortfall is the long tail of
analysis-specific column names that appear once or twice. Those are not left undefined by
oversight; a hand-written sentence for each is how this document's counts went stale twice.

| item | file |
|---|---|
| Table S1, panel and ENCODE accessions | `results/tables/supplementary_table_s1.csv` |
| Table S2, achieved match quality | `results/tables/match_quality.csv`, `match_quality_per_dataset.csv` |
| Table S3, fold integrity of retained scores, all three arms | `results/tables/fold_integrity.csv`, `fold_integrity_per_dataset.csv` |
| Table S4, standalone model AUROCs by arm | `results/tables/standalone_auroc.csv`, `standalone_auroc_per_dataset.csv` |
| Table S5, region asymmetry and the region-matched arm | `results/tables/region_asymmetry.csv`, `region_asymmetry_per_dataset.csv` |
| Table S6, ENCODE peak thresholds across the panel | `results/tables/peak_thresholds.csv`, `peak_thresholds_per_dataset.csv` |
| Table S7, design-effect components | `results/tables/design_effect.csv`, `design_effect_per_dataset.csv` |
| Table S8, positive-set overlap between arms | `results/tables/positive_set_overlap.csv` |
| Per-dataset results, three protocols by three models | `results/tables/three_arm_models_per_dataset.csv` |
| Supplementary figures S1 to S10 | `results/figures/f0,f1,f2,f3,f4,f5,f6,f7,f8,f13*` (every PDF in `results/figures/` that is not one of the main figures listed above) |
| Legends for the MAIN display items | typeset in `manuscript/paper.pdf` as figure and table captions |

| Supplementary figures S1 to S10, with legends | `manuscript/supplementary.pdf` |
| Supplementary Table S1 | `manuscript/supplementary_table_s1.csv` |

**The supplement is one document.** `manuscript/supplementary.tex` builds
`manuscript/supplementary.pdf`: ten figures numbered S1 to S10, each with a legend giving the
panels, the sample size, the uncertainty definition and the committed table its values come
from, plus a mapping from S-number to the repository build name and a note on Table S1's one
blank row. It is built by `manuscript/build.sh` alongside the paper, so it cannot go stale
against a regenerated figure.

Earlier releases shipped those ten as loose PDFs named `f0` to `f8` and `f13`, with no
S-numbering and no legends, and said so accurately. An audit was right that accuracy about an
unpublishable supplement is not a supplement. `tests/unit/test_supplement.py` now requires the
mapping to match the tree, the S-numbers to be contiguous, every legend to state its sample
size, uncertainty and source, and the main and supplementary figures to partition the built set
exactly, because five figures were orphaned once before when a section was cut.

## Reproducing the numbers

```
python scripts/verify.py --local results/tables
```

1156 numeric assertions against `config/golden.yaml`, and the number of assertions that ran is
itself asserted, so a check cannot silently skip. A clean `git clone` of this repository passes
all of them; that is the property worth checking, rather than that they pass in a working copy.

Two assertions are stronger than regression gates. `scripts/recompute.py` rebuilds 285 published
AUROCs from committed per-example scores to a maximum absolute difference of 3.3e-16, and
`scripts/k_sweep.py` rebuilds the headline contrast from raw sequence to 1.2e-06. Per-window
out-of-fold scores are committed for all three model classes and all three protocols
(`data/evidence/`), so every cell of the model-class comparison is recomputable here.

The full pipeline from raw ENCODE files needs a genome, cloud credentials and roughly 60 US
dollars of compute at current prices, and is not required to check any published value.
`run.sh` documents it and `docs/COST.md` breaks the figure down.

## Still outstanding

1. **Zenodo DOI.** The manuscript currently gives the GitHub URL only. `docs/ZENODO.md` has the
   procedure, and it is the single source of truth for it: this file used to give a different
   instruction, `docs/ZENODO.md` gave two different ones in two places, and
   `data-availability.tex` gave a fourth, so four documents described three workflows.
   When the DOIs exist, uncomment the sentence at the end of
   `manuscript/sections/data-availability.tex` and insert **both**: the **version** DOI and its
   tag, which identify the exact snapshot the reported numbers came from and are what a
   reproducibility citation needs, and the **concept** DOI beside it for discovery. Then run
   `cd manuscript && ./build.sh`. **Two places, not one**: the same two DOIs also go in
   **Code availability** in `paper.tex`, which currently gives only a moving GitHub URL. This
   item previously said Data availability was the only manuscript edit required, which
   contradicted `docs/ZENODO.md` and was wrong; both sections carry a commented slot.
2. **A journal submission after the preprint** will need the abstract cut to about 250 words
   and a main/supplement split. Neither is required by bioRxiv, and both depend on which
   journal, so neither is done.
3. **Per-commit AI co-author trailers stop after `202953d`; disclosure moves to the
   manuscript.** Most commits up to that point carry a `Co-Authored-By: Claude` trailer; the
   count is in `results/tables/history_scan.csv`, generated, because a number typed here is a
   number that is wrong by the next commit. That is accurate but it is the wrong granularity: it makes the commit list, which is the first thing
   a visitor to the repository sees, read as though no person was involved, and it undersells
   the author's own contribution. From `202953d` the trailer is no longer added. **This is a
   change of granularity and not of disclosure**: the historical trailers are left in place, no
   history is rewritten, and the manuscript's Use of AI tools section remains the disclosure of
   record. Recorded here so the change is visible rather than looking like something that
   quietly stopped.
4. **The AI-use disclosure needs an author decision before submission, and is deliberately not
   rewritten here.** The paragraph in `paper.tex` names one system and describes coding and
   drafting assistance. Multiple external audits produced by other generative systems were run
   against this repository during preparation, and their findings materially shaped it: they
   found defects that changed code, changed released tables and changed statements in the
   manuscript. Whether that is disclosable, and which systems to name, is the author's call
   under the target venue's policy, and it is an inventory only the author can make. A comment
   above the paragraph in `paper.tex` says the same thing so the decision is made rather than
   defaulted into. The sentence to check hardest is "no reference was generated by a language
   model", because it is a personal warranty: keep it only if it can be verified by hand
   against the bibliography.
5. **The historical billing-account ID needs an explicit decision, recorded on the release
   commit.** `SECURITY.md` states the finding and the three options.
   `scripts/history_scan.py --check` is a gate now, so the finding cannot fade out of the
   record between here and the tag, and it will still be reported on the commit that is
   archived. Deferred to the release phase on purpose; not dropped.

The reference list is verified, the figures are final, and the verifier and test suite pass on
a clean clone. This list previously read "Nothing else", which was wrong when it was written:
an external review then found eight stale counts in this file alone. Counts stated here are now
checked by `scripts/release_consistency.py`, which derives each one from the built artefact and
fails on a mismatch, so the way this section goes stale next will not be a number.

## What the last review round changed

Four independent reviews were run against the manuscript and the repository IN THAT ROUND.
They are enumerated in `docs/AI_USE_INVENTORY.md`, which is the one place that list lives; no
running total is restated here, because a count that rises with every audit goes stale in the
document that quotes it and did. None of the four broke a headline claim; all four returned major
revision on presentation and disclosure. The
substantive change is a new Results subsection: the bias-aware protocol matches fold only,
while both composition-matched protocols also match transcript region, so region alone
separates its classes at a median AUROC of 0.748 against exactly 0.5000 in the other two. That
asymmetry was undisclosed. Rebuilding the arm with the donor draw stratified on region lowers its composition baseline
from 0.8248 to 0.8052 and its contribution from +0.0122 to +0.0092, so 46% of its baseline
excess over the GC arm is region mix. The arm still carries the highest baseline and the lowest
contribution of the three, so the ordering is not a region artefact, and the span widens from
5.42 to 7.20.

The other changes worth naming: the Methods stated the paired-variance argument backwards; the
design effect of 1.35 was an unsourced constant and is now measured at 1.15, so the published
figure is conservative; the ENCODE peak files turn out to be pre-thresholded, so the Limitations
conceded a flaw the study does not have; `\citet{demler2012}` on DeLong for nested models is now
cited and answered; and the title no longer generalises the inference the paper exists to
refute. Twenty-five new assertions gate the new evidence, including the bias-aware arm's fold
integrity, which the manuscript had asserted was clean without ever measuring it.
