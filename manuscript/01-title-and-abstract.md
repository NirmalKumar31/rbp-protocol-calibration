# Title, abstract, and front matter

**Target venue:** *NAR Genomics and Bioinformatics* (Methods / benchmarking). Abstract must be a
single paragraph of at most 200 words with no citations, no figure references and no
non-standard abbreviations.

---

## Title

> **Harder negative sets lower apparent AUROC and raise measured model contribution: a
> three-protocol calibration across 94 ENCODE eCLIP datasets**

**Why this title and not the earlier ones.** Two previous framings were rejected in review.
"There is no protocol-independent measure of what a sequence model contributes" is a universal
negative, which invites a reader to go find a counterexample, and three protocols cannot support
a universal quantifier. "What a sequence model adds over composition is set by the benchmark's
headroom" names a mechanism the paper's own results partly retire: the baseline-to-gain gradient
is a property of composition-matched negatives and is absent for the bias-aware arm, so
"headroom" is not the general mechanism. The surviving title states the measured reversal, which
is the paper's most citable and least contestable fact, and it makes the negative set the agent.

## Running title

> Negative-set choice calibrates measured model contribution

## Author

Nirmal Kumar Thirupallikrishnan, Northeastern University.
Correspondence: thirupallikrishnan.n@northeastern.edu

## Abstract (198 words)

> Sequence models for RNA-binding proteins are evaluated against negative windows, and how those
> windows are chosen is usually treated as an implementation detail. Across 94 paired ENCODE
> eCLIP datasets we hold the model, the positives, the chromosome-blocked folds and the estimator
> fixed and vary only the negative-set protocol, measuring the nested contribution of a 4-mer
> model over a 19-feature mono- and dinucleotide composition baseline refit under each protocol.
> Replacing GC-matched with dinucleotide-matched negatives lowers apparent AUROC by 0.110 in 94
> of 94 datasets while raising the measured contribution from 0.027 to 0.066. A bias-aware
> protocol drawing negatives from other proteins' binding sites is the easiest of the three
> discriminations and yields the least, 0.012. Apparent difficulty and measured contribution
> therefore move in opposite directions, spanning 5.4-fold, and that span is about three times
> the gap between a k-mer model and a convolutional network on the same datasets. The
> composition baseline each protocol leaves behind accounts for most of the variation. No
> rescaling supplies a
> transportable remedy: an exponent that equalises our three protocols exists, but the exponent
> fitted to an independent benchmark differs from it 2.4-fold. We recommend reporting the
> composition-only AUROC obtained under the same protocol beside every headline AUROC, and not
> comparing contributions measured under different protocols.

## Keywords

RNA-binding proteins; eCLIP; benchmarking; negative-set construction; sequence composition;
incremental value; AUROC

## What the abstract deliberately does not say, and why

1. **No decomposition into "compression" and "protocol effect".** That quantity is not
   identified for any of the three model classes (Results, "Why we do not decompose the range"),
   and the raw contrast needs no transplant, no link function and no transportability
   assumption.
2. **No claim that the two-number report makes contributions comparable.** It does not: the
   recommendation fails its own pre-registered test on an independent benchmark. The abstract
   claims only that the baseline should be reported and that contributions should not be
   compared across protocols, both of which the evidence supports.
3. **No "invariant across model classes".** Model class is a small but detectable component of
   the multiplier, and SpliceBERT's is significantly lower.
4. **"Measured contribution" everywhere, never "true contribution", "hidden", "concealed",
   "inflated" or "proper".** Those words presuppose a protocol-free quantity the paper argues
   does not travel. The register is calibration, not bias.
5. **The magnitudes are stated as what they are, increments over a mono- and dinucleotide
   baseline.** The Discussion states that for a 4-mer most of that increment is one further
   order of composition. Burying that would be the single easiest thing for a referee to find.
