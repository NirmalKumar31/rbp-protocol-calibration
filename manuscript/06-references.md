# References

**Status.** Every entry below was verified against the publisher record or an equivalent
authoritative source on 2026-09-02. Entries that could not be verified were **removed rather
than estimated**, and the two claims that depended on them were rephrased so they no longer
need a citation. Nothing here is a reconstruction from memory.

1. ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the human genome.
   *Nature* 2012;**489**:57-74. doi:10.1038/nature11247

2. Van Nostrand EL, Pratt GA, Shishkin AA, *et al.* Robust transcriptome-wide discovery of
   RNA-binding protein binding sites with enhanced CLIP (eCLIP). *Nature Methods*
   2016;**13**:508-514. doi:10.1038/nmeth.3810

3. Van Nostrand EL, Freese P, Pratt GA, *et al.* A large-scale binding and functional map of
   human RNA-binding proteins. *Nature* 2020;**583**:711-719. doi:10.1038/s41586-020-2077-3

4. Horlacher M, *et al.* A systematic benchmark of machine learning methods for protein-RNA
   interaction prediction. *Briefings in Bioinformatics* 2023;**24**(5):bbad307.
   doi:10.1093/bib/bbad307. PMID:37635383
   *Data:* Zenodo, doi:10.5281/zenodo.10600977 (md5 verified on download).

5. Frankish A, Carbonell-Sala S, Diekhans M, *et al.* GENCODE: reference annotation for the
   human and mouse genomes in 2023. *Nucleic Acids Research* 2023;**51**(D1):D942-D949.
   doi:10.1093/nar/gkac1071

6. Grimm DG, Azencott CA, Aicheler F, *et al.* The evaluation of tools used to predict the
   impact of missense variants is hindered by two types of circularity. *Human Mutation*
   2015;**36**(5):513-523. doi:10.1002/humu.22768

7. Whalen S, Schreiber J, Noble WS, Pollard KS. Navigating the pitfalls of applying machine
   learning in genomics. *Nature Reviews Genetics* 2022;**23**:169-181.
   doi:10.1038/s41576-021-00434-9

8. DeLong ER, DeLong DM, Clarke-Pearson DL. Comparing the areas under two or more correlated
   receiver operating characteristic curves: a nonparametric approach. *Biometrics*
   1988;**44**:837-845. doi:10.2307/2531595

9. Sun X, Xu W. Fast implementation of DeLong's algorithm for comparing the areas under
   correlated receiver operating characteristic curves. *IEEE Signal Processing Letters*
   2014;**21**(11):1389-1393. doi:10.1109/LSP.2014.2337313

10. Somers RH. A new asymmetric measure of association for ordinal variables. *American
    Sociological Review* 1962;**27**(6):799-811. doi:10.2307/2090408

11. Firth D. Bias reduction of maximum likelihood estimates. *Biometrika*
    1993;**80**(1):27-38. doi:10.1093/biomet/80.1.27

12. Mood C. Logistic regression: why we cannot do what we think we can do, and what we can do
    about it. *European Sociological Review* 2010;**26**(1):67-82. doi:10.1093/esr/jcp006

13. Alipanahi B, Delong A, Weirauch MT, Frey BJ. Predicting the sequence specificities of DNA-
    and RNA-binding proteins by deep learning. *Nature Biotechnology* 2015;**33**:831-838.
    doi:10.1038/nbt.3300

14. Maticzka D, Lange SJ, Costa F, Backofen R. GraphProt: modeling binding preferences of
    RNA-binding proteins. *Genome Biology* 2014;**15**:R17. doi:10.1186/gb-2014-15-1-r17

15. Chen K, Zhou Y, Ding M, *et al.* Self-supervised learning on millions of primary RNA
    sequences from 72 vertebrates improves sequence-based RNA splicing prediction. *Briefings
    in Bioinformatics* 2024;**25**(3):bbae163. doi:10.1093/bib/bbae163

## Software

16. Harris CR, Millman KJ, van der Walt SJ, *et al.* Array programming with NumPy. *Nature*
    2020;**585**:357-362. doi:10.1038/s41586-020-2649-2

17. Virtanen P, Gommers R, Oliphant TE, *et al.* SciPy 1.0: fundamental algorithms for
    scientific computing in Python. *Nature Methods* 2020;**17**:261-272.
    doi:10.1038/s41592-019-0686-2

18. Pedregosa F, Varoquaux G, Gramfort A, *et al.* Scikit-learn: machine learning in Python.
    *Journal of Machine Learning Research* 2011;**12**:2825-2830.

19. Paszke A, Gross S, Massa F, *et al.* PyTorch: an imperative style, high-performance deep
    learning library. In: *Advances in Neural Information Processing Systems 32*, 2019:8024-8035.

20. Hunter JD. Matplotlib: a 2D graphics environment. *Computing in Science and Engineering*
    2007;**9**(3):90-95. doi:10.1109/MCSE.2007.55

---

## Removed, and how the text was changed instead

Five entries in the previous draft could not be verified in this pass. Rather than publish a
citation I had not confirmed, each was removed and the sentence depending on it rewritten:

- **Two incremental-value references** (Pencina *et al.*; Janes, Longton and Pepe). The
  Introduction previously cited them for the history of the incremental-value estimand. It now
  states plainly that the quantity is the difference in out-of-fold AUROC between a baseline
  model and the same model plus the score, which needs no citation because the paper defines it
  operationally in Methods.
- **A cross-cell-type generalisation reference** (Schreiber *et al.*). The Introduction's point
  about the sensitivity of genomic evaluations to how the negative class is defined is carried
  by [6] and [7], both verified.
- **A binormal ROC textbook** (Pepe). Methods gives the transformation
  d' = sqrt(2) x Phi^-1(AUROC) explicitly, so it is self-contained.
- **A dinucleotide-preserving shuffle reference** (Altschul and Erickson). The Discussion's
  justification for stopping the composition baseline at order two is now made directly, that
  order two is what a dinucleotide-preserving shuffle would hold fixed, without attribution to
  a specific shuffle algorithm since no shuffle is used in this work.

## Notes on citation placement

- **[4] belongs in Introduction sentence two.** It owns the phenomenon and this paper is a
  calibration of it, not a rediscovery. Any framing that reads as a rediscovery is both unfair
  and, for a preprint that will be read by that group, needlessly antagonistic.
- **[6] frames the gap.** It is the closest prior art in structure: benchmark assembly producing
  circularity that inflates and reorders methods.
- **[8] and [9] are both required.** [8] is the estimator, [9] is the implementation actually
  run, and at roughly 900,000 pooled pairs per dataset the distinction is not cosmetic.
- **[15] is the headline evidence for the model-class section and the project had never cited
  it.** A pretrained model with no source is not reproducible.
- **[3] carries two jobs:** the eCLIP resource, and the source of the significance thresholds
  that Limitations notes we deliberately did not apply.
- **PyTorch version.** Currently unpinned in the GPU requirements file and recoverable only from
  the container base image. Pin it or state the version in Methods before posting.
