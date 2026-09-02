# References

**Status.** Entries 1 to 12 were verified against the publisher record or an equivalent
authoritative source on 2026-09-02, and the fields shown are as verified. Entries marked
**[CHECK]** are needed by claims in the text but were not verified in this pass; confirm volume,
page and year before typesetting rather than trusting them. NAR style is numbered sequential
citation, first three authors then *et al.*, full article titles.

## Verified

1. ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the human genome.
   *Nature* 2012;**489**:57-74. **[CHECK volume/pages]**

2. Van Nostrand EL, Freese P, Pratt GA, *et al.* A large-scale binding and functional map of
   human RNA-binding proteins. *Nature* 2020;**583**:711-719. doi:10.1038/s41586-020-2077-3

3. Horlacher M, Cantini G, Hesse J, *et al.* A systematic benchmark of machine learning methods
   for protein-RNA interaction prediction. *Briefings in Bioinformatics*
   2023;**24**(5):bbad307. PMID:37635383
   *Data:* Zenodo, doi:10.5281/zenodo.10600977 (md5 verified on download).
   **[CHECK author list beyond the first]**

4. Grimm DG, Azencott CA, Aicheler F, *et al.* The evaluation of tools used to predict the
   impact of missense variants is hindered by two types of circularity. *Human Mutation*
   2015;**36**(5):513-523. doi:10.1002/humu.22768

5. Whalen S, Schreiber J, Noble WS, Pollard KS. Navigating the pitfalls of applying machine
   learning in genomics. *Nature Reviews Genetics* 2022;**23**:169-181. **[CHECK volume/pages]**

6. Schreiber J, Singh R, Bilmes J, Noble WS. A pitfall for machine learning methods aiming to
   predict across cell types. *Genome Biology* 2020;**21**:282. **[CHECK]**

7. Pencina MJ, D'Agostino RB, Pencina KM, *et al.* Interpreting incremental value of markers
   added to risk prediction models. *American Journal of Epidemiology*
   2012;**176**(6):473-481. **[CHECK]**

8. Janes H, Longton G, Pepe MS. Accommodating covariates in receiver operating characteristic
   analysis. *The Stata Journal* 2009;**9**(1):17-39. **[CHECK]**

9. DeLong ER, DeLong DM, Clarke-Pearson DL. Comparing the areas under two or more correlated
   receiver operating characteristic curves: a nonparametric approach. *Biometrics*
   1988;**44**:837-845.

10. Sun X, Xu W. Fast implementation of DeLong's algorithm for comparing the areas under
    correlated receiver operating characteristic curves. *IEEE Signal Processing Letters*
    2014;**21**(11):1389-1393. doi:10.1109/LSP.2014.2337313

11. Alipanahi B, Delong A, Weirauch MT, Frey BJ. Predicting the sequence specificities of DNA-
    and RNA-binding proteins by deep learning. *Nature Biotechnology* 2015;**33**:831-838.
    doi:10.1038/nbt.3300

12. Chen K, Zhou Y, Ding M, *et al.* Self-supervised learning on millions of primary RNA
    sequences from 72 vertebrates improves sequence-based RNA splicing prediction. *Briefings
    in Bioinformatics* 2024;**25**(3):bbae163. doi:10.1093/bib/bbae163
    *Weights:* `multimolecule/splicebert`, MultiMolecule package. Cite the package alongside
    the paper, since the package is the actual weight source used here.

## Needed and not yet verified

13. Van Nostrand EL, Pratt GA, Shishkin AA, *et al.* Robust transcriptome-wide discovery of
    RNA-binding protein binding sites with enhanced CLIP (eCLIP). *Nature Methods*
    2016;**13**:508-514. **[CHECK]** Needed for the assay and for ENCODE's standard peak
    thresholds, which the Limitations reference.

14. Frankish A, Carbonell-Sala S, Diekhans M, *et al.* GENCODE: reference annotation for the
    human and mouse genomes and transcriptomes. *Nucleic Acids Research* 2023;**51**:D942-D949.
    **[CHECK]** Cite the release actually used, v45.

15. Firth D. Bias reduction of maximum likelihood estimates. *Biometrika*
    1993;**80**(1):27-38. **[CHECK]**

16. Mood C. Logistic regression: why we cannot do what we think we can do, and what we can do
    about it. *European Sociological Review* 2010;**26**(1):67-82. **[CHECK]** Load-bearing:
    this is the diagnosis of the log-odds sign reversal reported in the model-class section.

17. Somers RH. A new asymmetric measure of association for ordinal variables. *American
    Sociological Review* 1962;**27**(6):799-811. **[CHECK]** Cited for the affine control in
    the transform sweep.

18. Maticzka D, Lange SJ, Costa F, Backofen R. GraphProt: modeling binding preferences of
    RNA-binding proteins. *Genome Biology* 2014;**15**:R17. **[CHECK]** Cited in Limitations as
    an example of the sampler design being discussed.

19. Pepe MS. *The Statistical Evaluation of Medical Tests for Classification and Prediction.*
    Oxford University Press, 2003. **[CHECK]** Cited for the binormal model underlying
    d' = sqrt(2) x Phi^-1(AUROC).

20. Altschul SF, Erickson BW. Significance of nucleotide sequence alignments: a method for
    random sequence permutation that preserves dinucleotide and codon usage. *Molecular Biology
    and Evolution* 1985;**2**(6):526-538. **[CHECK]** Cited for the dinucleotide-preserving
    shuffle that motivates stopping the composition baseline at order two.

## Software

21. Harris CR, Millman KJ, van der Walt SJ, *et al.* Array programming with NumPy. *Nature*
    2020;**585**:357-362.

22. Virtanen P, Gommers R, Oliphant TE, *et al.* SciPy 1.0: fundamental algorithms for
    scientific computing in Python. *Nature Methods* 2020;**17**:261-272.

23. Pedregosa F, Varoquaux G, Gramfort A, *et al.* Scikit-learn: machine learning in Python.
    *Journal of Machine Learning Research* 2011;**12**:2825-2830.

24. Paszke A, Gross S, Massa F, *et al.* PyTorch: an imperative style, high-performance deep
    learning library. *Advances in Neural Information Processing Systems* 2019;**32**:8024-8035.

25. Hunter JD. Matplotlib: a 2D graphics environment. *Computing in Science and Engineering*
    2007;**9**(3):90-95.

**Version pinning.** PyTorch is currently unpinned in the GPU requirements file and its version
is recoverable only from the container base image. Pin it or state the version explicitly in
Methods before submission.

## Position statement to cite in the AI disclosure

26. Committee on Publication Ethics. *Authorship and AI tools: COPE position statement.* 2023.
    Available from publicationethics.org. NAR's instructions direct authors to this statement.

## Notes on citation placement

- **[3] belongs in Introduction sentence two.** It owns the phenomenon and this paper is a
  calibration of it, not a rediscovery. Any framing that reads as a rediscovery is both unfair
  and a desk-rejection risk.
- **[4] frames the gap.** It is the closest prior art in structure: benchmark assembly producing
  circularity that inflates and reorders methods. Its absence from the Introduction would be the
  first thing a referee notes.
- **[9] and [10] are both required.** [9] is the estimator, [10] is the implementation actually
  run, and the distinction matters at this sample size.
- **[12] is the headline evidence in the model-class section and was uncited anywhere in the
  project until now.** A pretrained model with no source is not reproducible.
- **[2] carries two jobs**: the eCLIP resource, and the source of the significance thresholds
  that the Limitations note we deliberately did not apply.
