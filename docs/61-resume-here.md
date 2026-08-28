# Resume here

Updated 2026-08-28. **210/210 verifier checks, 594 tests, 9 figures, 2 manuscript orphans.**
Branch `r1-scale-check-and-corrected-refit`, **not merged to `master`**. Total spend ~$38.03.

## THE DECISION (council rounds 8 and 9): PAPER A

"Negative-set matching changes what an RBP benchmark measures." R1 plus five controls, at
*NAR Genomics & Bioinformatics* (Methods).

**Paper B is dead and must not be attempted.** Each of its three "traps" has a published owner
and two of the remedies are better than mine: Pencina, D'Agostino, Pencina, Janssens &
Greenland 2012 *Am J Epidemiol* 176(6):473-481 (AUROC-scale compression); Karlson, Holm & Breen
2012 *Sociological Methodology* 42:286-313 (cross-fit coefficient comparison, with the
significance test mine lacks); Schuster, Twisk, ter Riet, Heymans & Rijnhart 2021 *BMC Med Res
Methodol* 21:136 (the calibrated attenuation null; verified in Europe PMC). Whalen, Schreiber,
Noble & Pollard 2022 *Nat Rev Genet* 23:169-181 owns the "bundle of genomics ML pitfalls" slot.
Cite them; do not compete with them.

## WHAT THE PAPER IS

| section | claim | figure |
|---|---|---|
| **R1** | nested contribution +0.0265 GC-matched vs +0.0662 dinuc-matched; contrast **+0.0397** [+0.0336, +0.0458], 88/94; apparent AUROC falls 0.1095 in 94/94 | f0, f1 |
| **R1b** | 21% of the contrast is AUROC-scale compression; **+0.0313** survives, 87/94; d' scale +0.1290. Log-odds reversal reported and diagnosed | f8 |
| **R1c** | strand artifact bounded at **−0.0036** [−0.0071, +0.0001] against a region-matched placebo; **90.6%** survives | f3 |
| **R1d** | contrast replicates across cell lines at **r = +0.909**, better than either arm alone; efficiency **×1.31**, same conclusion on 58% of the data | f4 |
| **R1e** | rebuilt from raw sequence at k=4, difference **1.2e-06**; positive at every k=3..6 | f6 |
| **R1f** | **+0.0635** CDS-dominant vs **+0.0316** intron-dominant, p=1.5e-05; mechanism confirmed (intronic sites more compositional, p=2.7e-08) | f7 |
| R2 | methods table only, citing Horlacher 2023 as replication | f2 |

CUT: R3, R4/R4b/R4c, R4d, R5, R6, the composition-share framing. Replaced by one table saying
what was rejected and why.

## WHAT IS LEFT, IN ORDER

1. **Convert the audit document into a manuscript.** `docs/60-the-paper.md` is an internal
   source of truth organised as "R1, R1b, R1c...", not a preprint. It needs Title, Abstract,
   Introduction, Methods, Results, Discussion, References. Every number is settled and gated;
   none of it is written in the register a reader expects. **A day of writing, no analysis.**
   Referee #1's suggested abstract sentence is in the round-9 record.
2. **README rewrite.** What an interviewer actually opens.
3. **Publish.** Billing-ID scrub, git remote, push. Outward-facing and irreversible; needs
   explicit approval. Merge this branch to `master` first.

## STANDING CONSTRAINTS

Ask before Modal or GCS, every time. GCP GPU quota is 0. No em-dashes in the paper or slides.
Keep the Mac awake and cool during long runs (`caffeinate -is`, `OMP_NUM_THREADS=1`, `nice`).

## HARD-WON LESSONS THAT SHOULD NOT BE RELEARNED

- **A duplicated number is only as good as the assertion that the copies match.** `verify.py`
  read `cost_of_matching.csv` and never the two rehearsal tables it is a join of; a permuted
  dinuc table passed 166/166 while breaking 88/94, 94/94 and the Wilcoxon p.
- **Never compare a subset rebuild against a full-panel published mean.** That panel-mixing
  error appeared three separate times.
- **A restriction is not a random drop.** Any "restrict and recompute" control needs a matched
  placebo, and the placebo needs stratifying on whatever the restriction is correlated with.
- **Do not build a non-collapsibility null by drawing the covariate from the label** (collider),
  **or from a normal when the real covariate is skewed** (anti-conservative by ~10 points).
- **Gates written `if value is not None` disappear instead of failing.** `min_domain_checks`
  catches the whole class; it has caught me three times.
- **Golden keys built with an f-string are invisible** to `test_golden_keys_are_read`.
- Every one of the 210 assertions checks a *value*; none checked *which model produced it*,
  which is how the paper called a 4-mer a "5-mer" in four places and passed 150/150.
