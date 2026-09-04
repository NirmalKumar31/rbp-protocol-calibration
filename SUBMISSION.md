# Submission package

Everything needed to post the preprint and to reproduce every number in it. This file is the
index; nothing here is a summary of the science, which is in `manuscript/paper.pdf`.

## What to upload

| bioRxiv asks for | file |
|---|---|
| Manuscript PDF | `manuscript/paper.pdf` (28 pages) |
| Abstract (paste into the form) | abstract of `manuscript/paper.tex`, 246 words, no markup |
| Supplementary tables | `results/tables/supplementary_table_s1.csv` and the per-dataset tables listed below |
| Source, if requested | `manuscript/` is self-contained: `paper.tex`, `sections/`, `figures/`, `build.sh` |

**Subject area:** Bioinformatics. **Licence:** CC BY 4.0, matching the licence on `results/`
and `data/evidence/`. **Type:** New Results.

## The manuscript

`manuscript/paper.tex` builds `paper.pdf` with two runs of `pdflatex`, or in one step:

```
cd manuscript && ./build.sh
```

`build.sh` copies the six main figures out of `results/figures/` so the directory is a
self-contained upload, and it **fails** rather than warns on an undefined reference or citation.

Structure follows the convention of the target literature (Horlacher *et al.* 2023 and Chen
*et al.* 2024, both *Briefings in Bioinformatics*): Abstract, Introduction, Materials and
Methods, Results with declarative subsection headings, Discussion with Limitations, Conclusion,
then the declarations and 20 references.

`manuscript/` contains only what the submission needs: `paper.tex`, `sections/`, `figures/`,
`build.sh` and the built `paper.pdf`. Drafting notes and the record of editorial decisions are
on the `working-notes` branch.

## Main display items

| item | file | section |
|---|---|---|
| Figure 1 | `results/figures/f10_three_protocols.pdf` | three protocols, and the reversal |
| Figure 2 | `results/figures/f11_scale_sweep.pdf` | eight reparameterisations |
| Figure 3 | `results/figures/f12_protocol_or_baseline.pdf` | baseline attribution |
| Figure 4 | `results/figures/f9_deep_contrast.pdf` | model classes |
| Figure 5 | `results/figures/f14_external_validation.pdf` | independent benchmark |
| Figure 6 | `results/figures/f15_recommendation.pdf` | the recommendation, in and out of sample |

Tables 1 to 6 are typeset in the manuscript from the committed tables below. All figure PDFs
embed TrueType fonts and rasters are 400 dpi.

## Supplementary material

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
| Supplementary figures S1 to S10 | `results/figures/f0,f1,f2,f3,f4,f5,f6,f7,f8,f13*` (every PDF in `results/figures/` that is not one of the six main figures) |
| Legends for every display item | typeset in `manuscript/paper.pdf` as figure and table captions |

## Reproducing the numbers

```
python scripts/verify.py --local results/tables
```

674 numeric assertions against `config/golden.yaml`, and the number of assertions that ran is
itself asserted, so a check cannot silently skip. A clean `git clone` of this repository passes
all of them; that is the property worth checking, rather than that they pass in a working copy.

Two assertions are stronger than regression gates. `scripts/recompute.py` rebuilds 285 published
AUROCs from committed per-example scores to a maximum absolute difference of 2.2e-16, and
`scripts/k_sweep.py` rebuilds the headline contrast from raw sequence to 1.2e-06. Per-window
out-of-fold scores are committed for all three model classes and all three protocols
(`data/evidence/`), so every cell of the model-class comparison is recomputable here.

The full pipeline from raw ENCODE files needs a genome, cloud credentials and roughly 50 US
dollars of compute, and is not required to check any published value. `run.sh` documents it.

## Still outstanding

1. **Zenodo DOI.** `manuscript/sections/data-availability.tex` and the Code availability section
   currently give the GitHub URL only. See `docs/ZENODO.md` for the procedure; the DOI has to
   replace the placeholder and the manuscript then needs one rebuild.
2. Nothing else. The reference list is verified, the figures are final, and the verifier and
   test suite pass on a clean clone.

## What the last review round changed

Four independent reviews were run against the manuscript and the repository. None broke a
headline claim; all four returned major revision on presentation and disclosure. The
substantive change is a new Results subsection: the bias-aware protocol matches fold only,
while both composition-matched protocols also match transcript region, so region alone
separates its classes at a median AUROC of 0.748 against exactly 0.5000 in the other two. That
asymmetry was undisclosed. Rebuilding the arm with region matched lowers its composition
baseline from 0.8248 to 0.8017 and its contribution from +0.0122 to +0.0062, and the arm still
carries the highest baseline and the lowest contribution of the three, so the ordering is not a
region artefact and the span widens rather than narrowing.

The other changes worth naming: the Methods stated the paired-variance argument backwards; the
design effect of 1.35 was an unsourced constant and is now measured at 1.15, so the published
figure is conservative; the ENCODE peak files turn out to be pre-thresholded, so the Limitations
conceded a flaw the study does not have; `\citet{demler2012}` on DeLong for nested models is now
cited and answered; and the title no longer generalises the inference the paper exists to
refute. Twenty-five new assertions gate the new evidence, including the bias-aware arm's fold
integrity, which the manuscript had asserted was clean without ever measuring it.
