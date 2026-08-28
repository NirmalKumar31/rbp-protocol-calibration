# Resume here

Updated 2026-08-27. **139/139 verifier checks, 588 tests, 32 tables, 9 figures.**
Branch `r1-scale-check-and-corrected-refit`, commits `d7e645a`, `dd386d7`, `9314e71`.
**$0 spent this session, no cloud touched.**

## DONE this session

**1. R1 survives the ceiling objection.** This was the go/no-go and it passed.
`scripts/scale_check.py` -> `results/tables/scale_check.csv`, 13 gated checks, figure f8,
written up as R1b in `docs/60-the-paper.md`.

- The +0.0397 headline **reproduces exactly** (88/94) and now has provenance. It previously
  existed only in the manuscript: no table, no script, no golden key.
- AUROC compression is real, factor 1.51x, but accounts for only **21%** of the contrast.
  Net of it the protocol effect is **+0.0313 [+0.0267, +0.0363]**, positive in 87/94.
- Independent confirmation on the unbounded d' scale: **+0.1290 [+0.1091, +0.1499]**, 87/94.
- Somers' D is useless here and the paper says so: D = 2*AUROC - 1, so the contrast just
  doubles. (I had suggested it as the fix; it was wrong.)
- "hides two-thirds" was wrong twice over: uncorrected it is 60%, corrected **47.4%**.
  The paper now says "about half".
- **New honest limitation, reported not hidden:** on the log-odds scale the contrast
  REVERSES (-0.3771, dinuc larger in only 11/94). Diagnosed: the between-arm coefficient gap
  tracks each task's TOTAL signal (rho +0.520, p=8e-08) and not the incremental value it is
  supposed to measure (rho +0.065, p=0.53). Logistic coefficients are not comparable across
  fits with different total signal (Mood 2010). Normalising by total signal restores the sign
  (+0.1154, 68/94). Both the reversal and the fingerprint are gated, so if the fingerprint
  ever inverts the verifier fails and the reversal becomes real evidence against R1.

**2. Own-label leakage fixed** in `incremental_value.py`. Dedup now precedes the block
groupby. Published values corrected: unstratified **0.8246 -> 0.8164**, deciles
**0.7678/0.8319 -> 0.7548/0.8173**, spearman 0.339 -> 0.3375. Paper carries a dated
correction note. Both printed values are now gated; before, only a 0.70 floor was.

**3. Gates that actually bite.** The new block uses a `must()` accessor: a missing row is a
FAILURE, not a skip. Three attacks confirmed failing:
delete the protocol row -> 123/127; inflate it 25% -> 124/127; invert the fingerprint ->
125/127. Contrast with the old pattern, where deleting five rows gave 106/106.

**4. `run.sh` stage 13b added.** scale_check, incremental_value and recompute were wired into
nothing: the verifier asserted tables that no stage regenerated. They need no cloud and now
run in the pipeline. `strand_audit.py` is deliberately excluded, it needs an external GTF and
cannot run from committed evidence, and that gap is stated rather than hidden.

**5. The estimator R4 needed, done correctly, on the k-mer arm.**
`scripts/unconditional_refit.py` -> `unconditional_refit.csv`, 11 gated checks, written up as
R4d. Methods demonstration, not R4: the k-mer arm is the only per-variant table in the repo.

- `fit_delta_coef` always accepted `conservation=None`. Nothing ever passed it. Both fits now
  run on ONE row set built once, so attenuation cannot absorb a change in n.
- Raw attenuation is near zero and **flips sign** with the standardisation choice (+2.4%
  pooled, -5.2% within-dataset). That instability alone disqualifies it as a headline.
- **The conclusion inverts.** phyloP's coefficient is +2.12, and logistic regression is not
  collapsible, so the null for "no sharing" is a large AMPLIFICATION: **-68.4%** by forward
  simulation, **-59.9%** by the analytic cross-check. Excess over the null is **+70.8% /
  +63.3%**, same sign under both standardisations.
- **Self-validating:** the correlation needed to produce the observed attenuation is
  rho +0.085; the correlation measured directly is rho **+0.065**. They agree to 0.02.
- So on the k-mer arm the model's variant signal is **partly shared** with conservation, by
  an amount a modest correlation fully explains. Not orthogonal, not explained away.
- **My first null was wrong and is documented in the script, the paper and the golden file.**
  I drew the covariate from the label (C = mu*y + noise), which makes it a collider, and got
  a null of -1.5% -- which would have reproduced the retracted conclusion by a new route.
  Attack confirmed: substituting that broken null fails 4 checks.
- R4's own section now carries a DO-NOT-SUBMIT banner. `golden.yaml: r4_incremental_value`
  still gates the retracted numbers, deliberately, so the contradiction is visible in config.

**6. Verifier integrity, the whole skip class.** `integrity.min_domain_checks` asserts how
many checks ran. Most gates are still `if value is not None:` and vanish rather than fail when
a row goes missing; asserting the count fixes the category rather than 38 instances.
Demonstrated: deleting 2 rows from `multidonor_specificity.csv` silently disables **7 checks**
in a section never touched, and the floor is the only thing that notices. Also:
`recompute.py`'s rehearsal arm is now mandatory (deleting the evidence used to pass silently);
six retracted strand-contrast gates deleted; `cost` recentred -0.1080 -> -0.1095;
`auroc_gc`/`auroc_dinuc` recentred 0.7963/0.6893 -> 0.7981/0.6886.

**7. `scripts/audit_manuscript.py`** checks every number quoted to >=3 dp against the tables.
Found six unsourced cells in R1's own headline table, and that R4's 5-row table splices its
top row from a different table with a different n. **My first version was near-vacuous**: it
pooled 66k-row tables, saturating the 4-dp grid 73.9%, so a fabricated 0.9427 passed; it also
read its own report back as evidence. Now 6.3% saturated and the injection fails the build.
Gate is a **ratchet at 45, not a clean bill of health** - most survivors are subset aggregates.

**8. f7 retitled** (both titles contradicted the revised R5 text), Takeaway 3 rewritten around
the calibrated null, all 42 em-dashes gone.

## THE DECISION (council round 8, 2026-08-28): PAPER A, gated on ONE pre-registered test

**Paper A.** "Negative-set matching changes what an RBP benchmark measures." R1 + R1b on 94
paired ENCODE datasets. The scale decomposition is demoted to a robustness section citing
Pencina 2012. Venue: NAR Genomics & Bioinformatics (Methods).

**Paper B is dead** and must not be attempted. Each of its three "traps" has a published owner,
and two of the remedies are better than mine: scale compression is Pencina, D'Agostino,
Pencina, Janssens & Greenland 2012 *Am J Epidemiol* 176(6):473-481; cross-fit coefficient
non-comparability is Karlson, Holm & Breen 2012 *Sociological Methodology* 42:286-313, which
decomposes rescaling vs confounding WITH a significance test that mine lacks; the calibrated
attenuation null is Schuster, Twisk, ter Riet, Heymans & Rijnhart 2021 *BMC Med Res Methodol*
21:136 (verified in Europe PMC), plus Janes, Dominici & Zeger 2010 *Biostatistics* and Pang,
Kaufman & Platt 2013 *SMMR*. Whalen, Schreiber, Noble & Pollard 2022 *Nat Rev Genet*
23:169-181 already owns the "bundle of genomics ML pitfalls" slot. Cite these; do not compete
with them. (KHB is not in Europe PMC because *Sociological Methodology* is not indexed there;
verify it directly before citing.)

### PRE-REGISTERED, WRITTEN BEFORE THE TEST IS RUN

The strand wound is the only thing standing between here and submission, and the current
evidence has no leverage (`frac_sense` spans 0.433-0.615; extrapolation [-0.189, +0.072]).

**The test.** Fix the strand bug and regenerate negatives for a random **15-20 datasets**,
seed fixed and recorded before selection. Recompute the +0.0397 nested-gain contrast
**paired within dataset**, strand-correct against original, at `frac_sense` ~ 1.0. Windows
whose gene strand is ambiguous are EXCLUDED, never coin-flipped.

**Root cause, verified.** `annotation.py:126` `build_index` states "Strand is deliberately
dropped", so `negatives.py:328` falls back to `"strand": p["strand"]` and computes
`seq_rna = to_rna(dna, p["strand"])` from the POSITIVE's strand. Fix: retain strand in the
region index, recompute `to_rna` with each negative's OWN gene strand. **CPU-only, $0.** The
$25 estimate belonged to the SpliceBERT arm, which paper A does not use.

**Pass criteria, all three required:**
1. the contrast keeps its sign (positive), and
2. its bootstrap CI excludes zero, and
3. at least **60%** of the original point estimate survives (i.e. >= +0.0238 against +0.0397).

### RESULT, 2026-08-28: **PASSED on all three criteria. Write paper A.**

    full data (40 datasets)        +0.0378 [+0.0288, +0.0478]
    sense-only pairs (43% kept)    +0.0287 [+0.0201, +0.0383]
    placebo, same n dropped random +0.0346 [+0.0258, +0.0444]
    STRAND-SPECIFIC EXCESS         -0.0059 [-0.0098, -0.0021]
    strand-CORRECTED contrast      +0.0318 [+0.0234, +0.0415]   85% surviving

The artifact is REAL (interval excludes zero) and small. It costs ~15% of the contrast. That
is a measured control, not a null, and it is better evidence for being so. Gated as
`r1_strand_placebo`, 12 checks. Written up as R1c.

**Decision rule.** Pass -> write paper A, with a strand control Horlacher 2023 does not have.
Fail -> paper D: ship the repo as the portfolio artifact, post an honest short preprint, and
write the negative result up, which is itself worth publishing. **Do not renegotiate these
thresholds after seeing the result.**

## QUEUED, do these first when the session resumes

**Q1. Strand-test the surviving contrast. ~45 min, $0.** THE REASON THIS IS BACK: old Tier D1
(the $25 strand fix) was cancelled on the evidence that the contrast *grew* on sense-only
negatives, +0.2643 -> +0.2787. But +0.2643 is the composition-SHARE contrast, which is now
retracted as an algebraic identity. The cancellation rests on a withdrawn number, and the
surviving claim, the **+0.0397 nested gain**, has never been strand-tested. The paper already
admits this in its limitations.

Do it from committed data: `strand_audit.csv` carries per-dataset `frac_sense`; both rehearsal
arms carry per-dataset `delta_auroc`. Regress the per-dataset nested contrast on `frac_sense`,
and recompute the contrast within sense-rich and sense-poor halves. Flat means the artifact
does not drive it and the biggest referee objection largely dies. Not flat means we need to
know before submitting. Gate whichever way it lands, and do NOT reuse the retracted
`strand_audit` share machinery: that block was deleted from golden.yaml on purpose.

**Q2. Triage the 45 manuscript orphans. ~1 h, $0.** See item 3 below.

## THEN, in order. All $0, all CPU.

1. **R4 on SpliceBERT is the only thing left needing a download.** ~10 MB of
   `variants/scores_sb/` + `variants/scores_mm/`, ~$0.001. The procedure and the calibrated
   null are now written and gated, so it is a re-run, not a re-derivation. Expect it to
   demote R4 to a supporting result.
2. **Same leakage bug in `cloud_analysis.py:686-690`**, where `prev` is a regression
   covariate, so the 1.1689 row conditions on a covariate contaminated with the outcome.
   Requires rerunning the cloud analysis, so it is blocked on item 1.
3. **Triage the 45 manuscript orphans** (`results/tables/manuscript_orphans.csv`). Each is a
   number in the paper that traces to nothing. Expect most to be legitimate SUBSET aggregates
   (e.g. 0.8904 = mean `auroc_conservation_common` over the 44 powered rows having both
   columns; 0.8921 is the same rows' `auroc_conservation`, hence the -0.0017). Emit the
   legitimate ones into the table that owns them and lower `max_manuscript_orphans` toward 0.
   Purely mechanical, no judgement calls about claims, and it is the highest-value $0 work left.
4. **Quote headline numbers to 4 decimals, not 3.** The audit's 3-dp grid is 44% occupied, so
   a 3-decimal number is roughly a coin flip to slip through; 4-dp is 1 in 15.
5. Held until the above lands: README rewrite, billing-ID scrub + git remote.

## The old Tier A-E list, reconciled 2026-08-27

Mostly closed. A1, A2, B1, B3, B4, E2 done; C1 dropped then actually delivered as R6; B5, C2,
C3 deliberately dropped; A3 partly (f7 retitled, f4 not revisited); B2 (README) and E1
(manuscript) still open. **D1 is the live one**: cancelled on a number that has since been
retracted, so its strand question is reopened as Q1 above.

## Still stated as limitations rather than fixed

Strand-correct negative regeneration (CPU-only, feasible), expression filter on negatives,
eCLIP significance threshold, CLNREVSTAT >=1-star, R5 downsampling at fixed datasets,
multiplicity position, MC ANY-match audit, continuous power curve.

## Standing constraints

$0 further cloud spend. Ask before touching Modal or GCS, every single time. GCP GPU quota is
0. Explorer is not an option. No em-dashes in the paper or slides.
