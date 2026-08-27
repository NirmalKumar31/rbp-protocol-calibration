# Resume here

Updated 2026-08-27, after the R1 scale check and the corrected refit. **138/138 verifier
checks, 587 tests, 31 tables, 9 figures.** Uncommitted: last commit is still `60a5ded`.
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

## NEXT, in order. All $0, all CPU.

1. **R4 on SpliceBERT is the only thing left needing a download.** ~10 MB of
   `variants/scores_sb/` + `variants/scores_mm/`, ~$0.001. The procedure and the calibrated
   null are now written and gated, so it is a re-run, not a re-derivation. Expect it to
   demote R4 to a supporting result.
2. **Same leakage bug in `cloud_analysis.py:686-690`**, where `prev` is a regression
   covariate, so the 1.1689 row conditions on a covariate contaminated with the outcome.
   Requires rerunning the cloud analysis, so it is blocked on item 1.
3. **Remaining gate fixes:** delete the retracted strand contrast gate
   (`golden.yaml:396-398`, `verify.py:608-616`); add `min_verify_checks` to `main()`; convert
   the remaining `if v(k) is not None` skips to failures; make the rehearsal block mandatory
   in `recompute.py:134`; `cost` gate is centred on -0.1080 while the paper claims -0.1095
   (documented as spanning two panels, so check before changing).
4. **Manuscript:** retitle `figures.py` f7a/f7b, both still contradict the revised text.
   Rewrite Takeaway #3. Scrub em-dashes from the paper (new R1b text is already clean).
5. Held until the above lands: README rewrite, billing-ID scrub + git remote.

## Still stated as limitations rather than fixed

Strand-correct negative regeneration (CPU-only, feasible), expression filter on negatives,
eCLIP significance threshold, CLNREVSTAT >=1-star, R5 downsampling at fixed datasets,
multiplicity position, MC ANY-match audit, continuous power curve.

## Standing constraints

$0 further cloud spend. Ask before touching Modal or GCS, every single time. GCP GPU quota is
0. Explorer is not an option. No em-dashes in the paper or slides.
