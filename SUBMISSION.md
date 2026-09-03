# Submission package

Everything needed to post the preprint and to reproduce every number in it. This file is the
index; nothing here is a summary of the science, which is in `manuscript/paper.pdf`.

## What to upload

| bioRxiv asks for | file |
|---|---|
| Manuscript PDF | `manuscript/paper.pdf` (22 pages) |
| Abstract (paste into the form) | abstract of `manuscript/paper.tex`, 212 words, no markup |
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
| Table S3, fold integrity of retained scores | `results/tables/fold_integrity.csv`, `fold_integrity_per_dataset.csv` |
| Table S9, positive-set overlap between arms | `results/tables/positive_set_overlap.csv` |
| Per-dataset results, three protocols by three models | `results/tables/three_arm_models_per_dataset.csv` |
| Supplementary figures S1 to S8 | `results/figures/f0,f1,f2,f3,f4,f6,f13*` |
| Legends for every display item | typeset in `manuscript/paper.pdf` as figure and table captions |

## Reproducing the numbers

```
python scripts/verify.py --local results/tables
```

649 numeric assertions against `config/golden.yaml`, and the number of assertions that ran is
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
