"""Stage 13. Assert the reproduction actually reproduced, and fail loudly if not.

WHY THIS IS THE LAST STAGE AND NOT AN AFTERTHOUGHT. A pipeline that runs to completion and
quietly produces different science is worse than one that crashes, because nobody diffs a
plausible table. Reproducibility you do not check is not reproducibility; it is a hope with
a Dockerfile.

Every expectation lives in config/golden.yaml with an explicit tolerance, chosen to absorb
what legitimately varies (panel size, BLAS thread order, bootstrap seed, GPU vs CPU
inference at 1e-4) and nothing more.

    python scripts/verify.py                    # read tables from GCS
    python scripts/verify.py --local results/tables
    python scripts/verify.py --strict            # warnings become failures

Exit 0 means the science reproduced. Exit 1 means it did not, and the report says which
claim broke.
"""

import argparse
import io
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

from rbp.utils import cloud as cloudcfg  # noqa: E402

checks = []


def record(ok, claim, got, want, note=""):
    checks.append((ok, claim, got, want, note))
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {claim:44} got {got:<22} want {want}" + (f"  {note}" if note else ""),
          flush=True)


def near(claim, got, spec):
    """Point value within tolerance."""
    if got is None or (isinstance(got, float) and np.isnan(got)):
        return record(False, claim, "MISSING", f"{spec['value']} +/- {spec['tol']}")
    ok = abs(got - spec["value"]) <= spec["tol"]
    record(ok, claim, f"{got:.4f}", f"{spec['value']:.4f} +/- {spec['tol']}")


def at_least(claim, got, floor):
    ok = got is not None and got >= floor
    record(ok, claim, f"{got:.4f}" if got is not None else "MISSING", f">= {floor}")


def at_most(claim, got, ceil):
    ok = got is not None and got <= ceil
    record(ok, claim, f"{got:.3g}" if got is not None else "MISSING", f"<= {ceil:g}")


class Tables:
    """Result tables, from GCS or a local directory, whichever the run produced."""

    def __init__(self, local=None):
        self.local = Path(local) if local else None
        self.bucket = None if local else cloudcfg.bucket()

    def get(self, name):
        # Separator from the extension. Two of the artefacts are manifests written as TSV
        # (donor_tasks.tsv, variant_tasks.tsv) and reading them with the comma default yields a
        # single-column frame whose every attribute access raises AttributeError -- which is
        # how the pair-identity check below failed on its first run with a confusing error
        # rather than a missing file.
        sep = "\t" if name.endswith((".tsv", ".tsv.gz")) else ","
        if self.local:
            p = self.local / name
            return pd.read_csv(p, sep=sep) if p.exists() else None
        b = self.bucket.blob(f"results/tables/{name}")
        if not b.exists():
            return None
        return pd.read_csv(io.BytesIO(b.download_as_bytes()), sep=sep)


def verify_r1(T, g):
    print("\nR1  the negative-set protocol effect")
    d = T.get("cost_of_matching.csv")
    if d is None:
        return record(False, "cost_of_matching.csv present", "MISSING", "the table")
    spec = g["r1_cost_of_matching"]
    near("mean AUROC, GC arm", d.auroc_gc.mean(), spec["auroc_gc"])
    near("mean AUROC, dinuc arm", d.auroc_dn.mean(), spec["auroc_dinuc"])
    near("cost of proper matching", d.cost.mean(), spec["cost"])

    # THE NESTED GAIN, NOT THE DIFFERENCE OF TWO STANDALONE AUROCs.
    #
    # This computed (auroc - composition_auroc): how much better the sequence model is than a
    # composition model fitted separately. That is not the claim. The claim is that the
    # sequence score adds information OVER AND ABOVE composition, which is a nested comparison
    # -- composition alone against composition plus the score -- and the rehearsal already
    # computes it as delta_auroc = with_score_auroc - composition_auroc, with a bootstrap CI,
    # a p-value and a per-dataset `helps` flag attached.
    #
    # The two disagree materially: naive gives +0.0154 -> +0.0607 (ratio 3.94), nested gives
    # +0.0265 -> +0.0662 (ratio 2.50). Both support the direction, but the verifier was
    # certifying 3.94x while the manuscript's framing quotes the nested figure. A gate that
    # blesses a number nobody claims is worse than no gate, because it reads as confirmation.
    gain_gc = d.delta_auroc_gc.mean()
    gain_dn = d.delta_auroc_dn.mean()
    near("nested gain over composition, GC", gain_gc, spec["gain_gc"])
    near("nested gain over composition, dinuc", gain_dn, spec["gain_dinuc"])
    # THE THESIS. If this does not hold the paper's central claim is false, regardless of
    # how close every other number lands.
    at_least("gain RATIO (the paper's thesis)", gain_dn / gain_gc if gain_gc else None,
             spec["gain_ratio_min"])
    # The gain must be significant, not merely bigger. `helps` is the rehearsal's own verdict
    # on the nested comparison for that dataset.
    if "helps_dn" in d.columns:
        at_least("nested gain significant, dinuc arm", float(d.helps_dn.mean()),
                 spec["gain_significant_on_min_fraction"])
    near("fraction of datasets falling", float((d.cost < 0).mean()),
         spec["fraction_datasets_falling"])
    from scipy.stats import wilcoxon
    at_most("paired Wilcoxon p", wilcoxon(d.auroc_gc, d.auroc_dn)[1], spec["wilcoxon_p_max"])


def verify_r2(T, g):
    print("\nR2  four models on identical splits")
    # NOT `a or b`: pandas raises ValueError on DataFrame truthiness, so the fallback
    # would crash the verifier rather than fall through.
    d = T.get("matched_four_models.csv")
    if d is None:
        d = T.get("matched95_four_models.csv")
    if d is None:
        return record(False, "four-model table present", "MISSING", "the table")
    spec = g["r2_four_models"]
    col = {"composition": "composition_auroc", "kmer": "kmer_auroc",
           "cnn": "cnn", "splicebert": "splicebert"}
    means = {}
    for k, c in col.items():
        if c not in d:
            record(False, f"column {c}", "MISSING", "present")
            continue
        means[k] = d[c].mean()
        near(f"mean AUROC, {k}", means[k], spec[k])
    order = spec["required_order"]
    vals = [means.get(k) for k in order]
    ok = all(v is not None for v in vals) and all(a < b for a, b in zip(vals, vals[1:]))
    record(ok, "model ordering holds", " < ".join(order) if ok else str(vals),
           "strictly increasing")
    if "splicebert" in means:
        frac = float((d.splicebert > d.composition_auroc).mean())
        near("SpliceBERT beats composition, fraction", frac,
             spec["splicebert_beats_composition_fraction"])


def verify_r3(T, g):
    print("\nR3  positional concentration")
    d = T.get("locality_ism.csv")
    if d is None:
        return record(False, "locality_ism.csv present", "MISSING", "the table")
    spec = g["r3_locality"]
    near("k-mer Gini, median", d.kmer_gini.median(), spec["kmer_gini_median"])
    near("SpliceBERT Gini, median", d.sb_gini.median(), spec["splicebert_gini_median"])
    diff = d.sb_gini - d.kmer_gini
    near("median difference", diff.median(), spec["median_difference"])
    at_least("fraction SpliceBERT higher", float((diff > 0).mean()),
             spec["fraction_splicebert_higher_min"])

    # The strong form of the claim needs per-dataset uncertainty. If gini_sd is absent the
    # run used an older probe and cannot support "reversed on none".
    if {"kmer_gini_sd", "sb_gini_sd", "n_windows"} <= set(d.columns):
        se = np.sqrt((d.kmer_gini_sd ** 2 + d.sb_gini_sd ** 2) / d.n_windows)
        z = diff / se
        at_most("datasets significantly REVERSED", int((z < -1).sum()),
                spec["n_significantly_reversed_max"])
    else:
        record(False, "gini_sd carried through", "MISSING",
               "needed for the 'reversed on none' claim")
    from scipy.stats import wilcoxon
    at_most("paired Wilcoxon p", wilcoxon(d.sb_gini, d.kmer_gini)[1], spec["wilcoxon_p_max"])


def verify_r4_paired(T, g):
    """R4 as the PAIRED per-dataset comparison. This is the reported result.

    The pooled ladder is still checked by verify_r4 below, but it is no longer the claim:
    pooling inflated the matched arm from 0.755 to 0.829 through between-dataset scale.
    """
    print("\nR4  specificity, paired per dataset (the reported result)")
    per = T.get("variant_specificity.csv")
    paired = T.get("variant_ladder_paired.csv")
    coef = T.get("variant_coefficients.csv")
    ctrl = T.get("variant_specificity_controls.csv")
    if paired is None or per is None:
        return record(False, "paired specificity tables present", "MISSING", "the tables")
    spec = g["r4_paired_specificity"]

    hl = paired[paired.min_pathogenic == spec["min_pathogenic_headline"]]
    if not len(hl):
        return record(False, "headline stratum present", "MISSING",
                      f">={spec['min_pathogenic_headline']} pathogenic")
    r = hl.iloc[0]
    near("datasets in headline stratum", float(r.n_datasets), spec["n_datasets_headline"])
    near("AUROC, right-protein head", r.matched, spec["matched_auroc"])
    near("AUROC, wrong-protein head", r.mismatched, spec["mismatched_auroc"])
    near("AUROC, conservation", r.conservation, spec["conservation_auroc"])
    near("specificity gap", r.specificity_gap, spec["specificity_gap"])
    at_most("specificity paired p", float(r.p_specificity), spec["specificity_p_max"])
    at_least("right protein wins, fraction",
             float(r.matched_wins) / float(r.n_datasets), spec["matched_wins_min_fraction"])
    near("conservation still leads by", r.conservation_lead, spec["conservation_lead"])

    # The unflattering stratum. Asserted so it cannot be quietly dropped from the paper.
    allq = paired[paired.min_pathogenic == 0]
    if len(allq):
        near("all-datasets gap (shows nothing, by design)",
             allq.iloc[0].specificity_gap, spec["all_datasets_gap"])

    # A real effect grows with power; a pooling artefact need not.
    from scipy.stats import spearmanr
    rho = spearmanr(per.n_pathogenic, per.auroc_matched - per.auroc_mismatched).statistic
    at_least("gap grows with statistical power (rho)", float(rho), spec["gap_power_rho_min"])

    if coef is not None:
        w = coef[coef.standardisation == "within_dataset"].set_index("arm")
        near("coefficient, right protein (within dataset)",
             float(w.loc["matched", "coef"]), spec["coef_matched_within"])
        near("coefficient, wrong protein (within dataset)",
             float(w.loc["mismatched", "coef"]), spec["coef_mismatched_within"])
        if spec.get("coef_cis_must_separate"):
            sep = float(w.loc["matched", "ci_low"]) > float(w.loc["mismatched", "ci_high"])
            record(sep, "coefficient CIs separate (specificity established)",
                   f"matched low {w.loc['matched','ci_low']:.3f} vs mismatched high "
                   f"{w.loc['mismatched','ci_high']:.3f}", "no overlap")
        at_least("inference is clustered", float(w.loc["matched", "n_clusters"]),
                 spec["min_clusters"])

    # THE TRIVIAL POSITIONAL BASELINE. Asserted with the sign that embarrasses us, because a
    # gate that only checks flattering numbers is not a gate.
    if "block_prevalence" in hl and pd.notna(r.get("block_prevalence")):
        near("TRIVIAL positional baseline AUROC", float(r.block_prevalence),
             spec["block_prevalence_auroc"])
        at_most("model does NOT beat the trivial baseline",
                float(r.model_minus_prevalence), spec["model_minus_prevalence_max"])

    # The composition share, which had no interval until it was bootstrapped.
    rob = T.get("robustness.csv")
    if rob is not None:
        rr = rob.set_index("check")
        pairs = [("composition share, GC-matched", "composition_share_gc"),
                 ("composition share, dinuc-matched", "composition_share_dinuc"),
                 ("composition share, GC minus dinuc", "composition_share_drop")]
        for k, gk in pairs:
            if k in rr.index:
                near(f"R1 {k}", float(rr.loc[k, "value"]), spec[gk])
        for k in ("R1 effect vs log10(dataset size)",
                  "specificity gap vs size, powered stratum"):
            if k in rr.index:
                at_most(f"size effect bounded: {k[:34]}",
                        abs(float(rr.loc[k, "value"])), spec["size_rho_max"])

        # THE SHARE IS MODEL-DEPENDENT, and this is now the paper's headline. These keys sat
        # in golden.yaml unread for a day after being added -- Bug 29 (dead golden config)
        # committed alongside the fix it documents.
        for k, gk in (("composition share vs k-mer", "composition_share_vs_kmer"),
                      ("composition share vs CNN", "composition_share_vs_cnn"),
                      ("composition share vs SpliceBERT", "composition_share_vs_splicebert"),
                      ("composition share, k-mer minus SpliceBERT",
                       "composition_share_model_gap")):
            if k in rr.index:
                near(f"R1 {k}", float(rr.loc[k, "value"]), spec[gk])

        # THE IDENTITY, asserted instead of the CI. share_m = C/gain_m with C fixed across
        # models, so the ratio of two shares must equal the inverse ratio of two gains. A
        # "CI excludes zero" gate here would assert something true by construction; this
        # asserts the structure, and fails if the estimator is ever changed.
        four = T.get("matched_four_models.csv")
        if four is not None and {"kmer_auroc", "splicebert", "composition_auroc"} <= set(four):
            gk = four.kmer_auroc.mean() - 0.5
            gs = four.splicebert.mean() - 0.5
            c = four.composition_auroc.mean() - 0.5
            if gk > 0 and gs > 0:
                at_most("share ratio IS the inverse gain ratio (a rescaling, not a finding)",
                        abs((c / gk) / (c / gs) - gs / gk),
                        spec["share_ratio_equals_inverse_gain_ratio_max_diff"])

        # THE CROSS-CHECK. "composition share vs k-mer" and "composition share, dinuc-matched"
        # are the SAME quantity computed from two different tables (matched_four_models.csv and
        # cost_of_matching.csv). If they ever disagree, one of the tables is wrong and the
        # headline is built on it.
        a_k, b_k = "composition share vs k-mer", "composition share, dinuc-matched"
        if a_k in rr.index and b_k in rr.index:
            at_most("k-mer share agrees across two tables",
                    abs(float(rr.loc[a_k, "value"]) - float(rr.loc[b_k, "value"])),
                    spec["kmer_share_cross_check_max_diff"])

    # THE FOUR ATTACKS. Each was a reason the specificity result might not be real.
    at_ = T.get("variant_specificity_attacks.csv")
    if at_ is not None:
        a = at_.set_index("attack")
        def _v(k):
            return float(a.loc[k, "value"]) if k in a.index else None
        if _v("gap rises with power (spearman)") is not None:
            at_least("attack 1: gap rises with power (spearman)",
                     _v("gap rises with power (spearman)"),
                     spec["attack_gap_power_spearman_min"])
        if _v("wrong-protein floor is flat across thresholds") is not None:
            at_most("attack 1: wrong-protein floor stays flat",
                    _v("wrong-protein floor is flat across thresholds"),
                    spec["attack_floor_range_max"])
        k2 = "specificity survives the trivial rule as a covariate"
        if _v(k2) is not None:
            record(_v(k2) == 1.0 or not spec["attack_survives_trivial_rule"],
                   "attack 2: survives the trivial rule as covariate",
                   str(a.loc[k2, "note"])[:60], "CIs separate")
        if _v("permutation null is centred at zero") is not None:
            at_most("attack 3: permutation null at zero",
                    abs(_v("permutation null is centred at zero")),
                    spec["attack_null_gap_abs_max"])
        for k, gk in (("trivial rule at 100 kb", "trivial_rule_100kb"),
                      ("trivial rule at 1000 kb", "trivial_rule_1mb"),
                      ("trivial rule at 10000 kb", "trivial_rule_10mb")):
            if _v(k) is not None:
                near(f"attack 4: {k}", _v(k), spec[gk])

    # The wrong-protein floor must NOT be explained by donor similarity, or the control is
    # contaminated and the whole result goes with it.
    if ctrl is not None:
        c = ctrl.set_index("check")
        if "floor vs donor cell line" in c.index:
            at_least("floor independent of donor cell line (p)",
                     float(c.loc["floor vs donor cell line", "p"]),
                     spec["floor_cell_line_p_min"])
        if "floor vs donor's own strength" in c.index:
            at_most("floor barely tracks donor strength (rho)",
                    abs(float(c.loc["floor vs donor's own strength", "value"])),
                    spec["floor_donor_strength_rho_max"])


def verify_r4(T, g):
    print("\nR4  the POOLED ladder (retained for comparison, NOT the reported result)")
    d = T.get("variant_ladder.csv")
    if d is None:
        return record(False, "variant_ladder.csv present", "MISSING",
                      "the table (stage 11 aggregate)")
    spec = g["r4_variant_ladder"]
    row = {r.arm: r for r in d.itertuples()}
    got = {}
    for arm, key in [("kmer", "kmer_auroc"), ("mismatched", "mismatched_head_auroc"),
                     ("matched", "matched_head_auroc"),
                     ("conservation", "conservation_auroc")]:
        r = row.get(arm)
        got[key] = float(r.auroc) if r is not None else None
        near(f"AUROC, {arm}", got[key], spec[key])
    order = spec["required_order"]
    vals = [got.get(k) for k in order]
    ok = all(v is not None for v in vals) and all(a < b for a, b in zip(vals, vals[1:]))
    record(ok, "ladder is monotone", " < ".join(order) if ok else str(vals),
           "strictly increasing")

    for arm, key in [("matched", "matched_coef"), ("mismatched", "mismatched_coef"),
                     ("kmer", "kmer_coef")]:
        r = row.get(arm)
        near(f"clustered coefficient, {arm}",
             float(r.coef) if r is not None and hasattr(r, "coef") else None, spec[key])
    m, mm = row.get("matched"), row.get("mismatched")
    if m is not None and mm is not None:
        at_least("matched minus mismatched (specificity)",
                 float(m.coef) - float(mm.coef), spec["matched_minus_mismatched_min"])
        if hasattr(m, "n_clusters"):
            at_least("inference is clustered", int(m.n_clusters), spec["min_clusters"])


def verify_multidonor(T, g):
    """The multi-donor control, and the two panels that disagree.

    This replaced verify_donor_overlap, which gated the RETRACTED single-donor stratification.
    That verifier asserted a gap on "clean" donors and a CEILING on the floor-vs-overlap
    correlation -- a ceiling on the very association the claim needed to exist. It passed
    71/71 while the claim underneath it was confounded, which is the failure this function is
    written not to repeat.

    So the gates here are deliberately two-sided: the powered stratum must show the effect AND
    the full panel must show its absence. Asserting only the flattering panel is what went
    wrong last time.
    """
    print("\nmulti-donor wrong-protein control  (both panels, including the null one)")
    r = T.get("multidonor_specificity.csv")
    if r is None:
        return record(False, "multidonor_specificity.csv present", "MISSING", "the table")
    spec = g["r4_multidonor"]

    # THE ATTACK THIS BLOCKS, and it defeated every other gate in this function.
    #
    # A reviewer permuted the `donor` column WITHIN each target in multidonor_pairs.csv --
    # every gap, floor and target unchanged, only the donor-to-gap correspondence destroyed --
    # re-ran the analysis, and this verifier passed 86/86. Worse, the gate written specifically
    # to catch Bug 32 ("gap is NOT driven by donor quality") is satisfied OPTIMALLY by noise,
    # because scrambling guarantees rho ~ 0. A gate asserting the ABSENCE of an association
    # cannot distinguish "confound removed" from "labels meaningless"; it needs a paired
    # identity check, which is this.
    pairs = T.get("multidonor_pairs.csv")
    tasks = T.get("donor_tasks.tsv")
    if pairs is not None and tasks is not None:
        got = set(zip(pairs.dataset, pairs.donor))
        want = set(zip(tasks.target, tasks.donor))
        record(got <= want and len(got) > 0,
               "every (target, donor) pair is one the draw actually assigned",
               f"{len(got - want)} unassigned", "0")

    # AND THE PAIR-SET CHECK ABOVE IS NOT ENOUGH, which is the point worth remembering.
    # Permuting donors WITHIN a target preserves the pair set exactly, so the subset test
    # passes on scrambled labels -- I wrote it, ran the attack, and it passed 103/103. The only
    # thing that distinguishes "this donor produced this floor" from "these floors were dealt
    # out at random" is recomputing the floor from that donor's own per-variant scores. So it
    # is recomputed, for a sample, from data/evidence/scores_md/.
    md = ROOT / "data" / "evidence" / "scores_md"
    if pairs is not None and md.exists():
        from sklearn.metrics import roc_auc_score
        worst, n = 0.0, 0
        for _, row in pairs.head(40).iterrows():
            tp, tc = row.dataset.split(":")
            dp, dc = row.donor.split(":")
            f = md / f"{tc}_{tp}__{dc}_{dp}.csv"
            if not f.exists():
                continue
            s = pd.read_csv(f).dropna(subset=["delta"])
            if s.label.nunique() < 2:
                continue
            worst = max(worst, abs(roc_auc_score(s.label, s.delta.abs()) - float(row.auroc_floor)))
            n += 1
        at_least("donor floors recomputed from per-variant scores", n, 30)
        at_most("max |recomputed - reported| wrong-protein floor", worst, 1.0e-09)

    s = r.set_index(["panel", "estimator"])

    def cell(panel, est, col="value"):
        try:
            return float(s.loc[(panel, est), col])
        except KeyError:
            return None

    ALL, POW = "all usable", "powered (n_path>=20)"

    n_pairs = cell(ALL, "mean gap (unadjusted)", "n_pairs")
    if n_pairs is not None:
        near("multi-donor (target, donor) pairs", n_pairs, spec["n_pairs"])
        near("targets, all usable", cell(ALL, "mean gap (unadjusted)", "n_targets"),
             spec["n_targets_all"])
        near("targets, powered", cell(POW, "mean gap (unadjusted)", "n_targets"),
             spec["n_targets_powered"])

    # THE GATE THAT WOULD HAVE CAUGHT THE ORIGINAL ERROR. In the single-donor arm this was
    # -0.533: the gap tracked how much weaker the donor's model was. It must stay near zero,
    # or the control is measuring model capacity again.
    rho = cell(POW, "gap vs donor advantage (spearman)")
    if rho is not None:
        at_most("gap is NOT driven by donor quality (rho)", abs(rho),
                spec["gap_vs_donor_advantage_rho_max"])

    # powered stratum: the effect is there
    v = cell(POW, "mean gap (unadjusted)")
    if v is not None:
        near("powered stratum gap", v, spec["powered_gap"])
        at_least("powered wins, fraction",
                 cell(POW, "mean gap (unadjusted)", "wins")
                 / cell(POW, "mean gap (unadjusted)", "n_targets"),
                 spec["powered_wins_min_fraction"])
        at_most("powered paired p", cell(POW, "mean gap (unadjusted)", "p"),
                spec["powered_p_max"])
    v = cell(POW, "intercept | advantage+size")
    if v is not None:
        near("powered intercept at zero donor advantage", v, spec["powered_intercept"])
        if spec["powered_intercept_ci_excludes_zero"]:
            record(cell(POW, "intercept | advantage+size", "ci_low") > 0,
                   "powered intercept CI excludes zero",
                   cell(POW, "intercept | advantage+size", "ci_low"), "> 0")
    v = cell(POW, "donors STRONGER than target")
    if v is not None:
        near("gap where the DONOR is the better model", v, spec["stronger_donor_gap"])
        w, n = (cell(POW, "donors STRONGER than target", c) for c in ("wins", "n_targets"))
        if w and n:
            at_least("stronger-donor wins, fraction", w / n,
                     spec["stronger_donor_wins_min_fraction"])

    # THE MECHANISM behind the two panels disagreeing: a generic floor should not care how
    # many pathogenic variants a dataset has, and the protein's own head should.
    v = cell(POW, "floor vs power (spearman)")
    if v is not None:
        at_most("wrong-protein floor is FLAT in power", abs(v), spec["floor_vs_power_rho_max"])
    v = cell(POW, "matched arm vs power (spearman)")
    if v is not None:
        at_least("matched arm is STEEP in power", v, spec["matched_vs_power_rho_min"])

    # THE FULL PANEL, WHICH SHOWS NOTHING, asserted with its unflattering value
    v = cell(ALL, "mean gap (unadjusted)")
    if v is not None:
        near("all-panel gap (the null)", v, spec["all_panel_gap"])
    v = cell(ALL, "intercept | advantage+size")
    if v is not None:
        near("all-panel intercept (the null)", v, spec["all_panel_intercept"])
        record(not (cell(ALL, "intercept | advantage+size", "ci_low") > 0) if not spec["all_panel_intercept_ci_excludes_zero"] else True,
               "all-panel intercept CI INCLUDES zero, as reported",
               cell(ALL, "intercept | advantage+size", "ci_low"), "<= 0")
    # and the specification fragility, also asserted rather than hidden
    v = cell(POW, "intercept | advantage+size+power")
    if v is not None:
        near("powered intercept once power is adjusted", v,
             spec["powered_intercept_with_power"])
        record(not (cell(POW, "intercept | advantage+size+power", "ci_low") > 0) if not spec["powered_intercept_with_power_ci_excludes_zero"] else True,
               "power-adjusted intercept CI INCLUDES zero, as reported",
               cell(POW, "intercept | advantage+size+power", "ci_low"), "<= 0")


def verify_scale_check(T, g):
    """R1: the protocol effect must survive removal of AUROC-scale compression.

    THIS IS THE GATE ON THE PAPER'S PRIMARY CLAIM, and until now that claim had no gate at
    all. +0.0397 was printed in the manuscript, appeared in no committed table, and was
    produced by no script that could be found. It could not be reproduced and it could not
    fail. Every other headline in this project that reached that state turned out to be wrong.

    Missing rows are FAILURES here, not skips. The rest of this file guards its rows with
    `if v(k) is not None`, and a referee deleted five rows from a table and got 106/106. A
    check that a corrupted input can switch off is not a check.
    """
    print("\nR1  scale check  (is the protocol effect real, or the AUROC ceiling?)")
    d = T.get("scale_check.csv")
    if d is None:
        return record(False, "scale_check.csv present", "MISSING",
                      "run scripts/scale_check.py")
    spec = g["r1_scale_check"]
    q = d.set_index("check")

    def must(k, col="value"):
        """Read a row, or fail loudly. Returns None only after recording the failure."""
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, col])

    for k, gk in (
            ("k-mer size, both arms", "kmer_k"),
            ("mean AUROC, composition alone, GC arm", "mean_composition_gc"),
            ("mean AUROC, composition alone, dinuc arm", "mean_composition_dn"),
            ("mean AUROC, score alone, GC arm", "mean_score_gc"),
            ("mean AUROC, score alone, dinuc arm", "mean_score_dn"),
            ("mean AUROC, composition + score, GC arm", "mean_full_gc"),
            ("mean AUROC, composition + score, dinuc arm", "mean_full_dn"),
            ("nested contribution, GC arm", "nested_gc"),
            ("nested contribution, dinuc arm", "nested_dn"),
            ("CONTRAST, AUROC scale (published headline)", "contrast_auroc"),
            ("contrast attributable to SCALE alone", "contrast_scale_only"),
            ("CONTRAST, protocol effect net of scale", "contrast_protocol"),
            ("CONTRAST, protocol effect, REVERSE transplant", "contrast_protocol_reverse"),
            ("compression, REVERSE transplant", "contrast_scale_only_reverse"),
            ("protocol effect, logit link", "contrast_protocol_logit"),
            ("protocol effect, logit link REVERSE", "contrast_protocol_logit_reverse"),
            ("scale-null residual under d' exponent 0.5", "logodds_residual_p05"),
            ("scale-null residual under d' exponent 1.5", "logodds_residual_p15"),
            ("CONTRAST, d-prime scale (unbounded)", "contrast_dprime"),
            ("fraction of the contribution hidden by GC matching, corrected",
             "hidden_fraction_corrected"),
            ("CONTRAST, log-odds scale (REVERSES)", "contrast_logodds"),
            ("CONTRAST, log-odds normalised by total signal",
             "contrast_logodds_normalised"),
            ("log-odds gap predicted by SCALE alone", "logodds_scale_null"),
            ("log-odds residual after the scale null", "logodds_residual")):
        val = must(k)
        if val is not None:
            near(k, val, spec[gk])

    # THE CLAIM. The protocol effect must be positive with its interval clear of zero.
    lo = must("CONTRAST, protocol effect net of scale", "ci_low")
    if lo is not None and spec["protocol_ci_must_exclude_zero"]:
        record(lo > 0, "protocol effect survives with CI excluding zero", f"{lo:+.4f}", "> 0")

    # BOTH transplant directions must keep the sign, or the decomposition is a choice of
    # direction rather than a measurement.
    rlo = must("CONTRAST, protocol effect, REVERSE transplant", "ci_low")
    if rlo is not None and spec["reverse_ci_must_exclude_zero"]:
        record(rlo > 0, "protocol effect survives the UNFAVOURABLE transplant direction",
               f"{rlo:+.4f}", "> 0")

    sm = must("smallest protocol effect across links and directions")
    if sm is not None:
        at_least("protocol effect survives EVERY transplant choice", sm,
                 spec["min_protocol_effect_any_transplant"])

    # ...and compression must be a minority of it, or the headline IS the artefact.
    share = must("scale share of the published contrast")
    if share is not None:
        at_most("AUROC compression is a minority of the contrast", share,
                spec["max_scale_share"])

    # The log-odds reversal is disqualified only by this fingerprint. Assert both halves: the
    # coefficient gap must track total signal, and must NOT track incremental value. If that
    # ever flips, the reversal is real evidence against R1.
    rt = must("spearman(coef gap, TOTAL-signal gap)")
    if rt is not None:
        at_least("coef gap tracks TOTAL signal (why it is not comparable)", abs(rt),
                 spec["min_rho_coef_gap_total"])
    ri = must("spearman(coef gap, INCREMENTAL-value gap)")
    if ri is not None:
        at_most("coef gap does NOT track incremental value", abs(ri),
                spec["max_rho_coef_gap_increment"])


def verify_incremental_value(T, g):
    """R4: the model's variant signal net of conservation and position.

    This is the primary ClinVar claim, and it is gated on the same day it was
    promoted -- because the last thing promoted to a headline without a gate turned out to be
    an algebraic identity that had already been retracted once.

    The gate is on the attenuation FRACTION, not on a difference. Two coefficients drifting
    together would satisfy an absolute check and would not satisfy this one, which is the whole
    point: the claim is that conservation does not explain the model away, and that is a
    statement about a ratio.
    """
    print("\nR4  incremental value over conservation and position")
    d = T.get("incremental_value.csv")
    if d is None:
        return record(False, "incremental_value.csv present", "MISSING",
                      "run scripts/incremental_value.py")
    spec = g["r4_incremental_value"]
    q = d.set_index("check")

    def v(k):
        """A MISSING ROW IS A FAILURE. This was `return ... if k in q.index else None`, and
        every caller below guarded on `is not None`, so a referee deleted the five R4 rows
        from the table and the verifier reported 106/106 with no complaint. Rows are the
        evidence; their absence is the loudest possible signal, not the quietest."""
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, "value"])

    for k, gk in (("coef, matched, unconditional", "coef_matched_unconditional"),
                  ("coef, matched, controls = conservation only", "coef_matched_given_phylop"),
                  ("coef, matched, controls = plus trivial window rule",
                   "coef_matched_given_both"),
                  ("coef, mismatched, controls = conservation only",
                   "coef_mismatched_given_phylop")):
        if v(k) is not None:
            near(k[6:], v(k), spec[gk])

    # THE CLAIM.
    a = v("attenuation fraction after conditioning on phyloP")
    if a is not None:
        at_most("conservation does not explain the model away (attenuation)", abs(a),
                spec["max_attenuation_fraction"])
    s = v("arms separated under control")
    if s is not None and spec["arms_must_not_overlap_under_control"]:
        record(s == 1.0, "right- and wrong-protein CIs separated under control",
               "separated" if s == 1.0 else "OVERLAPPING", "separated")

    # And the second baseline must be a SECOND baseline, not conservation wearing a hat.
    m = v("positional rule, MIN AUROC within a phyloP decile")
    if m is not None:
        at_least("positional rule survives within every phyloP decile", m,
                 spec["positional_rule_min_auroc_within_phylop_decile"])
    # The printed values, which the floor above did not pin. Both moved when own-label leakage
    # was removed from the prevalence baseline, so they are asserted rather than trusted.
    u = v("positional rule, unstratified")
    if u is not None:
        near("positional rule, unstratified", u, spec["positional_rule_unstratified_auroc"])
    r = v("spearman(positional rule, phyloP)")
    if r is not None:
        near("spearman(positional rule, phyloP)", r, spec["positional_rule_phylop_spearman"])


def verify_strand_contrast(T, g):
    """Does the strand artifact drive R1's surviving contrast? Recorded, not cleared.

    The deleted strand block gated a retracted quantity behind a ceiling loose enough to pass a
    real confound. This one pins the measured values and asserts only what the data support:
    the halves agree, and the association is not the known size effect relabelled. The adverse
    extrapolation is gated too, so it cannot quietly disappear from the limitations.
    """
    print("\nR1  strand contrast  (does the artifact drive +0.0397?)")
    d = T.get("strand_contrast.csv")
    if d is None:
        return record(False, "strand_contrast.csv present", "MISSING",
                      "run scripts/strand_contrast.py")
    spec = g["r1_strand_contrast"]
    q = d.set_index("check")

    def must(k, col="value"):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, col])

    n = must("mean frac_sense on these datasets")
    if n is not None:
        near("mean frac_sense", n, spec["mean_frac_sense"])
        record(float(q.loc["mean frac_sense on these datasets", "n"])
               == spec["n_datasets"]["value"], "audited datasets in the test",
               q.loc["mean frac_sense on these datasets", "n"], spec["n_datasets"]["value"])

    for k, gk in (("mean contrast on these datasets", "contrast_on_subset"),
                  ("spearman(frac_sense, CONTRAST)", "rho_frac_sense_contrast"),
                  ("PARTIAL rho(frac_sense, contrast | size)", "partial_rho_given_size"),
                  ("difference between halves", "halves_difference"),
                  ("EXTRAPOLATED contrast at frac_sense = 1.0",
                   "extrapolated_at_full_sense")):
        v = must(k)
        if v is not None:
            near(k, v, spec[gk])

    # The association must not be the size effect wearing a costume.
    raw = must("spearman(frac_sense, CONTRAST)")
    par = must("PARTIAL rho(frac_sense, contrast | size)")
    if raw is not None and par is not None:
        at_most("association is not the size effect relabelled", abs(raw - par),
                spec["max_raw_partial_gap"])

    # THE ONE REASSURING RESULT. Halves agree, and the gap is small against the +0.0397.
    hd = must("difference between halves")
    lo = must("difference between halves", "ci_low")
    hi = must("difference between halves", "ci_high")
    if hd is not None:
        at_most("sense-rich and sense-poor halves agree", abs(hd),
                spec["max_abs_halves_difference"])
    if lo is not None and hi is not None and spec["halves_ci_must_include_zero"]:
        record(lo <= 0 <= hi, "halves difference straddles zero",
               f"[{lo:+.4f}, {hi:+.4f}]", "includes 0")

    # And the reason the correlation is weak, asserted so nobody reads p > 0.05 as a result.
    lev = must("frac_sense leverage available")
    if lev is not None:
        at_most("frac_sense leverage is small (test is weak BY DESIGN)", lev,
                spec["max_frac_sense_leverage"])


def verify_region(T, g):
    """Heterogeneity by binding region, gated together with the limitation on it."""
    print("\nR1f  region heterogeneity")
    d = T.get("region_heterogeneity.csv")
    if d is None:
        return record(False, "region_heterogeneity.csv present", "MISSING",
                      "run scripts/region_heterogeneity.py")
    spec = g["r1_region"]
    q = d.set_index("check")

    def must(k, col="value"):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, col])

    for k, gk in (("contrast, cds-dominant datasets", "contrast_cds"),
                  ("contrast, intron-dominant datasets", "contrast_intron"),
                  ("contrast, utr3-dominant datasets", "contrast_utr3"),
                  ("CDS minus intron", "cds_minus_intron"),
                  ("composition alone, intron-dominant", "composition_intron"),
                  ("composition alone, cds-dominant", "composition_cds")):
        v = must(k)
        if v is not None:
            near(k, v, spec[gk])

    lo = must("CDS minus intron", "ci_low")
    if lo is not None and spec["cds_intron_ci_must_exclude_zero"]:
        record(lo > 0, "CDS-intron difference excludes zero", f"{lo:+.4f}", "> 0")

    ci, cc = must("composition alone, intron-dominant"), must("composition alone, cds-dominant")
    if ci is not None and cc is not None and spec["intron_must_be_more_compositional"]:
        record(ci > cc, "MECHANISM: intronic sites are more compositional",
               f"{ci:.4f} vs {cc:.4f}", "intron higher")

    # THE LIMITATION. It must NOT survive adjustment, and that is asserted.
    pr = must("...partialling out total nested gain")
    if pr is not None:
        at_most("region acts THROUGH effect size, not independently", abs(pr),
                spec["max_partial_rho"])


def verify_deep_contrast(T, g):
    """R1g: the contrast is not a property of the model class -- it GROWS with capacity.

    R1's limitation used to be that every number came from an L2-penalised logistic on 4-mer
    counts. This gate holds a 7,089-parameter CNN and a 19.7M-parameter fine-tuned SpliceBERT
    to exactly the standard R1b is held to: the contrast positive with its interval clear of
    zero, and the protocol effect surviving EVERY member of the transplant family rather than
    the flattering one.

    Missing rows are FAILURES, per verify_scale_check. A gate a corrupted table can switch
    off is not a gate.
    """
    print("\nR1g deep-model contrast  (is R1 a property of bags of k-mers?)")
    d = T.get("deep_contrast.csv")
    if d is None:
        return record(False, "deep_contrast.csv present", "MISSING",
                      "run scripts/deep_model_contrast.py")
    spec = g["r1g_deep_contrast"]
    q = d.set_index(["model", "quantity"])

    def must(model, quantity, col="value"):
        if (model, quantity) not in q.index:
            record(False, f"row present: {model}/{quantity}", "MISSING", "the row")
            return None
        return float(q.loc[(model, quantity), col])

    contrasts = {}
    for model in ("kmer", "cnn", "splicebert"):
        sm = spec[model]
        for quantity, gk in (("nested_gc", "nested_gc"), ("nested_dn", "nested_dn"),
                             ("contrast_auroc", "contrast_auroc"),
                             ("protocol_effect_min", "protocol_min"),
                             ("protocol_effect_max", "protocol_max")):
            val = must(model, quantity)
            if val is not None:
                near(f"{model}: {quantity}", val, sm[gk])

        # The panel must be the whole panel, not whatever happened to finish.
        n = must(model, "contrast_auroc", "n")
        if n is not None:
            record(int(n) == spec["n_datasets"], f"{model}: datasets analysed",
                   int(n), spec["n_datasets"])

        pos = must(model, "contrast_positive_datasets")
        if pos is not None:
            # TWO-SIDED. Inflation passed everywhere while these were floors.
            near(f"{model}: datasets with a positive contrast", float(pos),
                 sm["positive_datasets"])

        lo = must(model, "contrast_auroc", "ci_low")
        contrasts[model] = must(model, "contrast_auroc")
        if lo is not None and spec["contrast_ci_must_exclude_zero"]:
            record(lo > 0, f"{model}: contrast interval excludes zero", f"{lo:+.4f}", "> 0")

        pmin = must(model, "protocol_effect_min")
        if pmin is not None and spec["protocol_must_survive_every_transplant"]:
            record(pmin > 0, f"{model}: protocol effect survives EVERY transplant choice",
                   f"{pmin:+.4f}", "> 0")

    # THE LADDER. The claim is that the contrast grows with model capacity, so the published
    # 4-mer number is the conservative end. If this fails the paper says "survives", not
    # "grows", and the wording must change rather than the gate.
    if spec["additive_ordering_holds"] and None not in contrasts.values():
        order = [contrasts["kmer"], contrasts["cnn"], contrasts["splicebert"]]
        record(order[0] < order[1] < order[2],
               "additive contrast is ordered kmer < cnn < splicebert (a fact about the additive "
               "scale, NOT a capacity claim -- that reading is withdrawn)",
               " < ".join(f"{v:+.4f}" for v in order), "kmer < cnn < splicebert")

    # The strong form: the deep model's interval must clear the k-mer's point estimate.
    sb_lo = must("splicebert", "contrast_auroc", "ci_low")
    if (sb_lo is not None and contrasts.get("kmer") is not None
            and spec["splicebert_ci_must_exceed_kmer_point"]):
        record(sb_lo > contrasts["kmer"],
               "SpliceBERT contrast interval clears the k-mer point estimate",
               f"{sb_lo:+.4f}", f"> {contrasts['kmer']:+.4f}")

    # ROWS. 1.0 in the GC arm is exact: it is scored from the file it trained on, so anything
    # less means a fold vanished.
    cg = must("-", "min_row_coverage_gc")
    if cg is not None:
        at_least("GC arm: every window scored", cg, spec["min_row_coverage_gc"])
    cd = must("-", "min_row_coverage_dn")
    if cd is not None:
        at_least("dinuc arm: window coverage after post-sweep drift", cd,
                 spec["min_row_coverage_dn"])

    # THE LADDER STEPS, PAIRED. Marginal intervals overlap between the k-mer and the CNN, so
    # the ordering is only a claim as a paired difference on the same 94 datasets.
    for a, b, gk in (("cnn", "kmer", "step_cnn_minus_kmer"),
                     ("splicebert", "cnn", "step_splicebert_minus_cnn"),
                     ("splicebert", "kmer", "step_splicebert_minus_kmer")):
        val = must("ladder", gk)
        if val is not None:
            near(f"ladder step: {a} - {b}", val, spec["ladder"][gk])
        lo = must("ladder", gk, "ci_low")
        if lo is not None and spec["ladder_steps_ci_must_exclude_zero"]:
            record(lo > 0, f"ladder step {a} - {b}: interval excludes zero",
                   f"{lo:+.4f}", "> 0")
    # The weakest rung is asserted at its own floor so it cannot quietly erode.
    weak = must("ladder", "step_cnn_minus_kmer_datasets")
    if weak is not None:
        near("CNN beats the k-mer contrast on this many datasets", float(weak),
             spec["step_datasets_cnn_minus_kmer"])

    # THE RATIO SCALE, WHERE THE LADDER REVERSES. R1b's rule applied to R1g's own headline:
    # the additive ladder is real, the multiplier ladder is not, and both must be reported.
    rspec = spec["ratio"]
    nboth = must("ratio", "datasets_positive_both_arms_all_models")
    if nboth is not None:
        record(int(nboth) == rspec["datasets_positive_both_arms"],
               "datasets with positive gains in both arms, all models",
               int(nboth), rspec["datasets_positive_both_arms"])
    mults = {}
    for model, gk in (("kmer", "multiplier_kmer"), ("cnn", "multiplier_cnn"),
                      ("splicebert", "multiplier_splicebert")):
        val = must("ratio", f"multiplier_{model}")
        mults[model] = val
        if val is not None:
            near(f"protocol multiplier, {model}", val, rspec[gk])
            at_least(f"{model}: multiplier is substantial on the ratio scale", val,
                     spec["min_multiplier_any_model"])
    for a, b, gk in (("splicebert", "kmer", "logstep_splicebert_minus_kmer"),
                     ("splicebert", "cnn", "logstep_splicebert_minus_cnn")):
        val = must("ratio", f"logstep_{a}_minus_{b}")
        if val is not None:
            near(f"log-ratio step: {a} - {b}", val, rspec[gk])
        hi = must("ratio", f"logstep_{a}_minus_{b}", "ci_high")
        if hi is not None and spec["ratio_ladder_must_reverse_for_splicebert"]:
            record(hi < 0, f"ratio ladder REVERSES for splicebert against {b} "
                           f"(interval clear of zero, wrong side)", f"{hi:+.4f}", "< 0")

    # THE SUMMARY MUST BE ARITHMETIC ON THE PER-DATASET TABLE, not an independent assertion.
    # Everything above reads deep_contrast.csv alone, so editing that one file would pass
    # every check. Recompute the means from the evidence.
    per = T.get("deep_contrast_per_dataset.csv")
    if per is None:
        record(False, "deep_contrast_per_dataset.csv present", "MISSING", "the evidence table")
    else:
        worst = 0.0
        for model in ("kmer", "cnn", "splicebert"):
            gain_gc = per[f"{model}_gain_gc"]
            gain_dn = per[f"{model}_gain_dn"]
            for quantity, series in (("nested_gc", gain_gc), ("nested_dn", gain_dn),
                                     ("contrast_auroc", gain_dn - gain_gc)):
                got = must(model, quantity)
                if got is not None:
                    worst = max(worst, abs(got - float(series.mean())))
            pos = must(model, "contrast_positive_datasets")
            if pos is not None:
                record(int(pos) == int((gain_dn - gain_gc > 0).sum()),
                       f"{model}: positive-dataset count matches the evidence",
                       int(pos), int((gain_dn - gain_gc > 0).sum()))
        at_most("summary means are arithmetic on the per-dataset table", worst,
                spec["max_summary_arithmetic_diff"])

    # --- ANTI-FORGERY -------------------------------------------------------------------
    if per is not None:
        # 1. The internal identity DeLong guarantees. Exact.
        worst_id = 0.0
        for model in ("kmer", "cnn", "splicebert"):
            for arm in ("gc", "dn"):
                worst_id = max(worst_id, float(
                    (per[f"{model}_gain_{arm}"]
                     - (per[f"{model}_full_{arm}"] - per[f"comp_{arm}"])).abs().max()))
        at_most("gain equals full minus composition, exactly, in every cell", worst_id,
                spec["max_gain_identity_diff"])

        # 2. The panel must be the R1 panel, by name.
        reh = T.get("rehearsal_binding_gc.csv")
        if reh is None:
            record(False, "rehearsal_binding_gc.csv present for the panel cross-check",
                   "MISSING", "the table")
        elif spec["panel_must_match_rehearsal"]:
            got, want = set(per.dataset), set(reh.dataset)
            record(got == want, "R1g's panel is exactly R1's panel, by dataset name",
                   f"{len(got & want)} shared, {len(got - want)} extra",
                   f"all {len(want)}")

        # 3. Coverage columns must EXIST, not default.
        if spec["coverage_columns_must_be_present"]:
            for arm in ("gc", "dn"):
                record(f"coverage_{arm}" in per.columns,
                       f"coverage_{arm} column is present in the evidence "
                       f"(absence must not read as perfect coverage)",
                       f"coverage_{arm}" in per.columns, True)

        # 4. THE ANCHOR. Recompute every cell's raw pooled AUROC from the committed
        # per-window scores. This is the check that makes the table non-forgeable.
        from sklearn.metrics import roc_auc_score
        roots = {"gc": ROOT / "data" / "evidence" / "scores_gc",
                 "dn": ROOT / "data" / "evidence" / "scores"}
        worst_raw, n_cells = 0.0, 0
        for r in per.itertuples():
            for model in ("cnn", "splicebert"):
                for arm in ("gc", "dn"):
                    base = roots[arm] / r.cell / r.protein / model
                    parts = []
                    for f in range(5):
                        fp = base / f"fold{f}" / "scores.tsv.gz"
                        if not fp.exists():
                            parts = None
                            break
                        parts.append(pd.read_csv(fp, sep="\t"))
                    if parts is None:
                        continue
                    sc = pd.concat(parts, ignore_index=True)
                    claimed_n = getattr(r, f"{model}_nrows_{arm}", None)
                    if claimed_n is None or int(claimed_n) != len(sc):
                        record(False, f"{r.dataset} {model} {arm}: score row count",
                               len(sc), claimed_n)
                        continue
                    got = roc_auc_score(sc.label.values, sc.score.values)
                    worst_raw = max(worst_raw,
                                    abs(got - float(getattr(r, f"{model}_raw_{arm}"))))
                    n_cells += 1
        at_most("raw pooled AUROC recomputed from committed per-window scores", worst_raw,
                spec["max_raw_auroc_diff"])
        at_least("cells anchored to committed per-window evidence", n_cells,
                 spec["min_anchored_cells"])

    # THE TIE TO R1. The 4-mer refitted here must reproduce the published contrast, or the
    # comparison against the deep models is measuring two different things.
    sc = T.get("scale_check.csv")
    if sc is None:
        record(False, "scale_check.csv present for the R1g cross-check", "MISSING", "the table")
    elif contrasts.get("kmer") is not None:
        pub = sc.set_index("check")
        key = "CONTRAST, AUROC scale (published headline)"
        if key not in pub.index:
            record(False, "published contrast row present", "MISSING", "the row")
        else:
            drift = abs(contrasts["kmer"] - float(pub.loc[key, "value"]))
            at_most("refitted 4-mer reproduces the published R1 contrast", drift,
                    spec["max_kmer_refit_drift"])


def verify_protocol_identification(T, g):
    """R1h: the protocol effect is not identified for ANY model, and the CNN survives most.

    The honest output is a specification span, not a point estimate. Gated so that a future
    refactor cannot pick the most favourable slope source and call it survival.
    """
    print("\nR1h protocol identification  (is the protocol effect an estimand?)")
    d = T.get("protocol_identification.csv")
    sl = T.get("protocol_baseline_slopes.csv")
    if d is None or sl is None:
        return record(False, "protocol_identification.csv + slopes present", "MISSING",
                      "run scripts/protocol_identification.py")
    spec = g["r1h_protocol_identification"]
    q = d.set_index("check")

    def must(k, col="value"):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, col])

    s = sl.set_index(["model", "arm"])
    for model, arm, gk in (("kmer", "gc", "slope_kmer_gc"),
                           ("splicebert", "gc", "slope_splicebert_gc")):
        if (model, arm) in s.index:
            near(f"baseline-invariance slope, {model} {arm}",
                 float(s.loc[(model, arm), "slope_net_of_mechanism"]), spec[gk])
        else:
            record(False, f"slope row present: {model}/{arm}", "MISSING", "the row")
    for model, arm, gk in (("kmer", "gc", "mechanical_slope_kmer_gc"),
                           ("splicebert", "dn", "mechanical_slope_splicebert_dn")):
        if (model, arm) in s.index:
            near(f"mechanical slope, {model} {arm}",
                 float(s.loc[(model, arm), "mechanical_slope"]), spec[gk])
    if ("kmer", "gc") in s.index:
        near("corr(full, composition) estimate errors, kmer gc",
             float(s.loc[("kmer", "gc"), "corr_full_comp"]), spec["corr_full_comp_kmer_gc"])
    # THE TELL. Setting the covariance to zero makes this constant within an arm.
    if spec["mechanical_slope_must_vary_across_models"]:
        for arm in ("gc", "dn"):
            vals = sl[sl.arm == arm].mechanical_slope.round(9).nunique()
            record(vals > 1, f"mechanical slope varies across models ({arm} arm)", vals,
                   "> 1 distinct value")

    for model, lo_k, hi_k in (("kmer", "link_family_min_kmer", "link_family_max_kmer"),
                              ("splicebert", "link_family_min_splicebert",
                               "link_family_max_splicebert")):
        lo, hi = must(f"{model}/link_family_min"), must(f"{model}/link_family_max")
        if lo is not None:
            near(f"link family minimum, {model}", lo, spec[lo_k])
        if hi is not None:
            near(f"link family maximum, {model}", hi, spec[hi_k])
        if lo is not None and spec["link_family_must_keep_sign"]:
            record(lo > 0, f"{model}: sign holds across every ROC-motivated link",
                   f"{lo:+.4f}", "> 0")

    for model, gk in (("kmer", "odds_forward_kmer"),
                      ("splicebert", "odds_forward_splicebert")):
        v = must(f"{model}/odds_forward")
        if v is not None:
            near(f"odds link (no ROC rationale), {model}", v, spec[gk])
            if spec["odds_link_must_reverse"]:
                record(v < 0, f"{model}: the odds link REVERSES the sign, and is excluded "
                              f"from the headline range for a stated reason",
                       f"{v:+.4f}", "< 0")

    for model, lo_k, hi_k, n_k in (
            ("kmer", "spec_span_min_kmer", "spec_span_max_kmer", "spec_excluding_zero_kmer"),
            ("splicebert", "spec_span_min_splicebert", "spec_span_max_splicebert",
             "spec_excluding_zero_splicebert")):
        lo, hi = must(f"{model}/spec_span_min"), must(f"{model}/spec_span_max")
        if lo is not None:
            near(f"specification span minimum, {model}", lo, spec[lo_k])
        if hi is not None:
            near(f"specification span maximum, {model}", hi, spec[hi_k])
        n = must(f"{model}/spec_excluding_zero")
        if n is not None:
            record(int(n) == spec[n_k], f"{model}: specifications excluding zero",
                   int(n), spec[n_k])
    tot = must("kmer/spec_total")
    if tot is not None:
        record(int(tot) == spec["spec_total"], "specifications per model", int(tot),
               spec["spec_total"])

    ncnn = must("cnn/spec_excluding_zero")
    if ncnn is not None:
        record(int(ncnn) == spec["spec_excluding_zero_cnn"],
               "cnn: specifications excluding zero", int(ncnn),
               spec["spec_excluding_zero_cnn"])
        if spec["cnn_must_survive_all_specifications"]:
            record(int(ncnn) == int(tot or 6),
                   "the CNN (the MIDDLE rung) is the only model surviving every "
                   "specification, so no capacity ordering holds", int(ncnn), int(tot or 6))
    nsb = must("splicebert/spec_excluding_zero")
    if nsb is not None and spec["splicebert_must_fail_majority_of_specifications"]:
        record(nsb * 2 < (tot or 6),
               "SpliceBERT's protocol effect fails a majority of specifications",
               f"{int(nsb)}/{int(tot or 6)}", "fails > half")


def verify_expression_control(T, g):
    """R1j: the expression confound is huge, balanced across arms, and costs 10% of the contrast.

    The balance assertion is the load-bearing one. A confound present equally in both arms
    inflates both absolute AUROCs and cannot manufacture a difference between them; a
    differential one can. If `arm_difference_must_straddle_zero` ever fails, the contrast is
    contaminated and this control stops being a defence.
    """
    print("\nR1j expression control  (are the negatives merely untranscribed?)")
    d = T.get("expression_control.csv")
    if d is None:
        return record(False, "expression_control.csv present", "MISSING",
                      "run scripts/expression_control.py")
    spec = g["r1j_expression_control"]
    q = d.set_index("check")

    def must(k, col="value"):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, col])

    n = q["n"].iloc[0] if "n" in q.columns else None
    if n is not None:
        record(int(n) == spec["n_datasets"], "datasets with both arms and RNA-seq",
               int(n), spec["n_datasets"])

    for k, gk in (
            ("negatives in an untranscribed locus, GC arm", "untranscribed_fraction_gc"),
            ("negatives in an untranscribed locus, dinuc arm", "untranscribed_fraction_dn"),
            ("contrast, full data", "contrast_full"),
            ("contrast, expressed-negative pairs", "contrast_expressed_only"),
            ("contrast, PLACEBO stratified on region", "contrast_placebo_stratified"),
            ("EXPRESSION-SPECIFIC EXCESS (stratified)", "expression_excess"),
            ("locus-mix component", "locus_mix"),
            ("expression-CORRECTED contrast", "contrast_corrected"),
            ("fraction of the contrast surviving", "surviving_fraction")):
        val = must(k)
        if val is not None:
            near(k, val, spec[gk])

    # THE BALANCE ARGUMENT, asserted both ways.
    lo = must("arm difference, untranscribed fraction", "ci_low")
    hi = must("arm difference, untranscribed fraction", "ci_high")
    if lo is not None and hi is not None and spec["arm_difference_must_straddle_zero"]:
        record(lo <= 0 <= hi, "the untranscribed fraction is BALANCED across arms, so the "
                              "confound cannot manufacture the contrast",
               f"[{lo:+.4f}, {hi:+.4f}]", "contains 0")
    diff = must("arm difference, untranscribed fraction")
    if diff is not None:
        at_most("absolute arm difference in untranscribed fraction", abs(diff),
                spec["max_abs_arm_difference"])

    # The artifact is real and small. Assert BOTH halves.
    ehi = must("EXPRESSION-SPECIFIC EXCESS (stratified)", "ci_high")
    if ehi is not None and spec["expression_excess_must_exclude_zero"]:
        record(ehi < 0, "the expression artifact is REAL (interval excludes zero), not noise",
               f"{ehi:+.4f}", "< 0")
    surv = must("fraction of the contrast surviving")
    if surv is not None:
        at_least("fraction of the contrast surviving the expression correction", surv,
                 spec["min_surviving_fraction"])
    noise = must("mean placebo seed noise (sd), GC arm")
    if noise is not None:
        at_most("placebo Monte Carlo noise stays below the effect measured", noise,
                spec["max_placebo_seed_noise"])


def verify_cluster_intervals(T, g):
    """R1i: the panel is 79 proteins, not 94 independent datasets."""
    print("\nR1i clustered intervals  (94 datasets are only 79 proteins)")
    d = T.get("cluster_intervals.csv")
    if d is None:
        return record(False, "cluster_intervals.csv present", "MISSING",
                      "run scripts/cluster_intervals.py")
    spec = g["r1i_cluster_intervals"]
    q = d.set_index("check")

    def must(k, col="value"):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, col])

    r = must("within_protein_correlation_kmer")
    if r is not None:
        near("within-protein correlation of the contrast, k-mer", r,
             spec["within_protein_correlation_kmer"])
    npairs = q.loc["within_protein_correlation_kmer", "n_pairs"] \
        if "within_protein_correlation_kmer" in q.index else None
    if npairs is not None and not pd.isna(npairs):
        record(int(npairs) == spec["n_doubled_proteins"],
               "proteins assayed in both cell lines", int(npairs),
               spec["n_doubled_proteins"])
    w = must("max_width_ratio")
    if w is not None:
        near("widest interval inflation under protein clustering", w,
             spec["max_width_ratio"])
        at_least("protein clustering actually widens the intervals", w,
                 spec["min_width_ratio"])
    for k, gk in (("R1_contrast", "r1_contrast_width_ratio"),
                  ("R1f_cds_minus_intron", "r1f_cds_width_ratio")):
        if k in q.index:
            near(f"{k}: interval inflation under protein clustering",
                 float(q.loc[k, "width_ratio"]), spec[gk])
            record(bool(q.loc[k, "excludes_zero_clustered"]),
                   f"{k} still excludes zero under protein clustering", True, True)
        else:
            record(False, f"row present: {k}", "MISSING", "the row")
    dup = must("R1c_duplicated_proteins")
    if dup is not None:
        record(int(dup) == spec["r1c_duplicated_proteins"],
               "R1c has no clustering to correct (its datasets are distinct proteins)",
               int(dup), spec["r1c_duplicated_proteins"])

    # The per-dataset counts, published and corrected. Both, so the change stays visible.
    for k, gk in (("helps_gc_published", "helps_gc_published"),
                  ("helps_gc_design_effect", "helps_gc_design_effect"),
                  ("helps_dn_published", "helps_dn_published"),
                  ("helps_dn_design_effect", "helps_dn_design_effect")):
        v_ = must(k)
        if v_ is not None:
            record(int(v_) == spec[gk], f"datasets where the score helps: {k}",
                   int(v_), spec[gk])

    flag = must("all_headlines_exclude_zero_clustered")
    if flag is not None and spec["all_headlines_exclude_zero_clustered"]:
        record(flag == 1.0,
               "every headline still excludes zero under protein clustering",
               bool(flag), True)


def verify_three_arm(T, g):
    """R1k: three protocols, a 5.4-fold range, and a falsified prediction kept on the record."""
    print("\nR1k three negative-set protocols  (is there a protocol-independent contribution?)")
    d = T.get("three_arm_contrast.csv")
    if d is None:
        return record(False, "three_arm_contrast.csv present", "MISSING",
                      "run scripts/three_arm_contrast.py")
    spec = g["r1k_three_arm"]
    q = d.set_index("check")

    def must(k, col="value"):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, col])

    comp, gain = {}, {}
    for arm, cg, gg in (("gc", "composition_gc", "gain_gc"),
                        ("dn", "composition_dn", "gain_dn"),
                        ("neg2", "composition_neg2", "gain_neg2")):
        c = must(f"composition alone, {arm} arm")
        v = must(f"nested contribution, {arm} arm")
        comp[arm], gain[arm] = c, v
        if c is not None:
            near(f"composition alone, {arm} arm", c, spec[cg])
        if v is not None:
            near(f"nested contribution, {arm} arm", v, spec[gg])

    for a, b, gk in (("dn", "gc", "contrast_dn_minus_gc"),
                     ("neg2", "gc", "contrast_neg2_minus_gc"),
                     ("neg2", "dn", "contrast_neg2_minus_dn")):
        v = must(f"CONTRAST, {a} minus {b}")
        if v is not None:
            near(f"contrast, {a} minus {b}", v, spec[gk])
    for a, b, gk in (("dn", "gc", "multiplier_dn_over_gc"),
                     ("neg2", "gc", "multiplier_neg2_over_gc")):
        v = must(f"protocol multiplier, {a} over {b}")
        if v is not None:
            near(f"protocol multiplier, {a} over {b}", v, spec[gk])

    # THE FALSIFIED PREDICTION, asserted in the direction the data actually went.
    if None not in comp.values() and spec["neg2_must_have_highest_composition"]:
        record(comp["neg2"] == max(comp.values()),
               "neg2 has the HIGHEST composition baseline (the pre-registered prediction "
               "said it would be the lowest)", f"{comp['neg2']:.4f}",
               f"> {max(comp['gc'], comp['dn']):.4f}")
    if None not in gain.values() and spec["neg2_must_have_lowest_gain"]:
        record(gain["neg2"] == min(gain.values()),
               "neg2 gives the LOWEST nested contribution of the three protocols",
               f"{gain['neg2']:+.4f}", f"< {min(gain['gc'], gain['dn']):+.4f}")
    if None not in gain.values():
        rng = max(gain.values()) / max(min(gain.values()), 1e-9)
        at_least("fold range in measured contribution across the three protocols", rng,
                 spec["min_fold_range_across_protocols"])
    cn = must("CONTRAST, neg2 minus gc")
    if cn is not None and spec["contrast_neg2_minus_gc_must_be_negative"]:
        record(cn < 0, "the field's own bias-aware protocol reveals LESS than GC matching, "
                       "so 'prefer harder negatives' is not supportable",
               f"{cn:+.4f}", "< 0")

    resid = []
    for src, tgt, gk in (("gc", "dn", "protocol_effect_gc_on_dn"),
                         ("gc", "neg2", "protocol_effect_gc_on_neg2"),
                         ("dn", "neg2", "protocol_effect_dn_on_neg2"),
                         ("neg2", "dn", "protocol_effect_neg2_on_dn")):
        v = must(f"protocol effect, {src} increment on {tgt} baseline")
        if v is not None:
            near(f"protocol effect, {src} -> {tgt}", v, spec[gk])
            resid.append(v)
    if resid and spec["transplant_residuals_must_change_sign"]:
        record(min(resid) < 0 < max(resid),
               "transplant residuals CHANGE SIGN with direction, so the decomposition does "
               "not identify a protocol effect",
               f"{min(resid):+.4f} to {max(resid):+.4f}", "spans zero")

    rho = must("spearman(composition baseline, nested gain), all arms pooled")
    if rho is not None:
        near("spearman(composition baseline, nested gain), pooled over arms", rho,
             spec["spearman_composition_gain"])


def verify_baseline_confounding(T, g):
    """R1l: the protocol/baseline distinction is not identifiable, and the paper says so."""
    print("\nR1l baseline confounding  (is 'protocol changes the contribution' trivial?)")
    d = T.get("baseline_confounding.csv")
    if d is None:
        return record(False, "baseline_confounding.csv present", "MISSING",
                      "run scripts/baseline_confounding.py")
    spec = g["r1l_baseline_confounding"]
    q = d.set_index("check")

    def must(k):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, "value"])

    for k, gk in (("R2, cubic in composition baseline alone", "r2_baseline_alone"),
                  ("R2, plus protocol dummies", "r2_with_protocol"),
                  ("F statistic, protocol beyond baseline", "f_protocol_beyond_baseline"),
                  ("common support width (AUROC)", "common_support_width"),
                  ("composition baseline median, gc arm", "baseline_median_gc"),
                  ("composition baseline median, dn arm", "baseline_median_dn"),
                  ("composition baseline median, neg2 arm", "baseline_median_neg2"),
                  ("fraction of gc->dn difference explained by compression",
                   "compression_share_gc_dn"),
                  ("fraction of gc->neg2 difference explained by compression",
                   "compression_share_gc_neg2"),
                  ("fraction of dn->neg2 difference explained by compression",
                   "compression_share_dn_neg2"),
                  ("R2, one constant increment on each cell's own baseline",
                   "r2_one_constant_increment"),
                  ("residual after the constant-increment model, dn arm", "residual_dn"),
                  ("residual after the constant-increment model, neg2 arm", "residual_neg2")):
        v = must(k)
        if v is not None:
            near(k, v, spec[gk])

    pv = must("p value, protocol beyond baseline")
    if pv is not None:
        at_least("the naive 'protocol beyond baseline' test stays UNCONVINCING, which is why "
                 "the overlap argument is needed", pv,
                 spec["min_p_protocol_beyond_baseline"])
    n = must("cells inside the common support")
    if n is not None:
        at_most("cells where two protocols share a composition baseline", int(n),
                spec["max_cells_in_common_support"])
    for k in ("fraction of gc->dn difference explained by compression",
              "fraction of gc->neg2 difference explained by compression",
              "fraction of dn->neg2 difference explained by compression"):
        v = must(k)
        if v is not None:
            at_most(f"compression is a minority of: {k.split(' difference')[0][12:]}", v,
                    spec["max_compression_share"])
    rd, rn = must("residual after the constant-increment model, dn arm"), \
        must("residual after the constant-increment model, neg2 arm")
    if rd is not None and rn is not None and spec["residuals_must_be_ordered"]:
        record(rd > 0 > rn,
               "residuals after the ceiling model are still ORDERED by protocol, so the arms "
               "differ by more than the arithmetic", f"dn {rd:+.4f}, neg2 {rn:+.4f}",
               "dn > 0 > neg2")

    # THE CLUSTERED SPEARMAN. Point estimates, and then the pattern of the intervals, which is
    # what the paper actually claims: the relation is present in both composition-matched arms
    # and not detectable in the bias-aware one. Gating the endpoints alone would pass a run in
    # which every interval had widened to include zero.
    # Golden keys spelled out and not built with an f-string: test_golden_keys_are_read.py
    # greps this file for each leaf name, and a constructed key is one nobody can find.
    ci = {}
    for arm, gspec in (("pooled", spec["spearman_clustered_pooled"]),
                       ("within gc arm", spec["spearman_clustered_gc"]),
                       ("within dn arm", spec["spearman_clustered_dn"]),
                       ("within neg2 arm", spec["spearman_clustered_neg2"])):
        k = f"spearman(baseline, gain) {arm}, protein-clustered CI"
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            continue
        near(k, float(q.loc[k, "value"]), gspec)
        ci[arm] = (float(q.loc[k, "ci_low"]), float(q.loc[k, "ci_high"]))
    if len(ci) == 4 and spec["spearman_clustered_signs_must_hold"]:
        excl = {a: hi < 0 for a, (lo, hi) in ci.items()}
        record(excl["pooled"] and excl["within gc arm"] and excl["within dn arm"]
               and not excl["within neg2 arm"],
               "clustered intervals: negative for pooled/gc/dn, spanning zero for neg2",
               ", ".join(f"{a.replace('within ', '').replace(' arm', '')} "
                         f"{'excl' if e else 'incl'}" for a, e in excl.items()),
               "pooled/gc/dn exclude 0, neg2 includes 0")


def verify_scale_sweep(T, g):
    """R1m: the search for a protocol-independent rescaling, and its failure."""
    print("\nR1m scale sweep  (does ANY rescaling remove the protocol dependence?)")
    d = T.get("scale_sweep.csv")
    if d is None:
        return record(False, "scale_sweep.csv present", "MISSING",
                      "run scripts/scale_sweep.py")
    spec = g["r1m_scale_sweep"]
    q = d.set_index("check")

    def must(k, col="value"):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, col])

    n = int(d.check.str.startswith("fold range,").sum())
    record(n == spec["n_transforms"], "transforms searched", n, spec["n_transforms"])
    for k, gk in (("fold range, raw AUROC gain", "fold_range_raw"),
                  ("fold range, d' increment (binormal)", "fold_range_dprime"),
                  ("fold range, logit increment", "fold_range_logit"),
                  ("fold range, headroom-normalised, g/(1-comp)", "fold_range_headroom"),
                  ("fold range, excess-normalised, g/(comp-0.5)", "fold_range_excess"),
                  ("minimum fold range over all transforms", "minimum_fold_range")):
        v = must(k)
        if v is not None:
            near(k, v, spec[gk])

    floor = must("minimum fold range over all transforms")
    if floor is not None:
        at_least("NO transform reaches protocol independence: the floor stays clear of 1",
                 floor, spec["min_floor_fold_range"])
    lo = must("fold range, headroom-normalised, g/(1-comp)", "ci_low")
    if lo is not None and spec["floor_ci_must_exclude_one"]:
        record(lo > 1.0, "the best coordinate's interval EXCLUDES 1, so it is a better "
                         "coordinate and not an invariant", f"{lo:.2f}", "> 1")
    raw = must("fold range, raw AUROC gain")
    som = must("fold range, Somers' D gain")
    if raw is not None and som is not None:
        at_most("Somers' D is a linear rescale so its fold range must EQUAL the raw one",
                abs(raw - som), spec["max_somers_raw_difference"])


def verify_protocol_or_baseline(T, g):
    """R1n: the protocol label adds ~1% once the baseline is known. It is the baseline."""
    print("\nR1n protocol or baseline  (does the protocol label carry information?)")
    d = T.get("protocol_or_baseline.csv")
    if d is None:
        return record(False, "protocol_or_baseline.csv present", "MISSING",
                      "run scripts/protocol_or_baseline.py")
    spec = g["r1n_protocol_or_baseline"]
    q = d.set_index("check")

    def must(k, col="value"):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, col])

    ip = must("incremental R2 of the protocol label, given the baseline")
    ib = must("incremental R2 of the baseline, given the protocol label")
    if ip is not None:
        near("incremental R2 of the protocol label", ip, spec["incremental_r2_protocol"])
        at_most("the protocol label stays uninformative given the baseline", ip,
                spec["max_incremental_r2_protocol"])
    if ib is not None:
        near("incremental R2 of the baseline", ib, spec["incremental_r2_baseline"])
    if ip and ib:
        at_least("the baseline explains far more than the protocol label does", ib / ip,
                 spec["min_baseline_over_protocol_ratio"])

    n = must("datasets where neg2 raises the composition baseline")
    if n is not None:
        record(abs((94 - int(n)) - spec["n_discordant_datasets"]["value"])
               <= spec["n_discordant_datasets"]["tol"],
               "datasets where neg2 LOWERS the baseline", 94 - int(n),
               spec["n_discordant_datasets"]["value"])
    con = must("neg2 minus gc gain, concordant datasets")
    dis = must("neg2 minus gc gain, discordant datasets")
    if con is not None:
        near("neg2 minus gc, concordant", con, spec["gain_diff_concordant"])
    if dis is not None:
        near("neg2 minus gc, discordant", dis, spec["gain_diff_discordant"])
    if con is not None and dis is not None and spec["discordant_must_reverse_sign"]:
        record(con < 0 < dis,
               "THE NATURAL EXPERIMENT: the sign REVERSES where neg2 lowers the baseline, so "
               "the baseline predicts and the protocol label does not",
               f"{con:+.4f} -> {dis:+.4f}", "negative -> positive")
    r = must("within-dataset spearman(delta baseline, delta gain)")
    if r is not None:
        near("within-dataset spearman(delta baseline, delta gain)", r,
             spec["within_dataset_spearman"])

    m = must("dn minus gc, matched on baseline")
    lo = must("dn minus gc, matched on baseline", "ci_low")
    hi = must("dn minus gc, matched on baseline", "ci_high")
    if m is not None:
        near("dn minus gc matched on baseline", m, spec["matched_dn_minus_gc"])
    if lo is not None and hi is not None and spec["matched_must_contain_zero"]:
        record(lo <= 0 <= hi,
               "the published contrast does NOT survive matching on the baseline",
               f"[{lo:+.4f}, {hi:+.4f}]", "contains 0")
    # THE CORRECTION. R1l's "no common support" was the three-way intersection; pairwise the
    # gc-vs-neg2 comparison is easy, and there a residual DOES exist.
    for label, key in (
            ("pairwise common support width, gc vs dn", spec["pairwise_support_gc_dn"]),
            ("pairwise common support width, gc vs neg2", spec["pairwise_support_gc_neg2"]),
            ("incremental R2 of protocol, gc vs dn", spec["incremental_r2_gc_dn"]),
            ("incremental R2 of protocol, gc vs neg2", spec["incremental_r2_gc_neg2"]),
            ("neg2 minus gc, matched on baseline", spec["matched_neg2_minus_gc"]),
            ("neg2 minus gc, within-dataset intercept at zero baseline shift",
             spec["intercept_neg2_minus_gc"]),
            ("within-arm spearman(baseline, gain), gc", spec["within_arm_spearman_gc"]),
            ("within-arm spearman(baseline, gain), neg2", spec["within_arm_spearman_neg2"])):
        v = must(label)
        if v is not None:
            near(label, v, key)
    mhi = must("neg2 minus gc, matched on baseline", "ci_high")
    ihi = must("neg2 minus gc, within-dataset intercept at zero baseline shift", "ci_high")
    if mhi is not None and ihi is not None and spec["neg2_residual_must_exclude_zero"]:
        record(mhi < 0 and ihi < 0,
               "a protocol-specific residual DOES exist for neg2 at matched baseline, under "
               "two independent designs", f"{mhi:+.4f}, {ihi:+.4f}", "both < 0")
    sg = must("within-arm spearman(baseline, gain), neg2")
    if sg is not None and spec["gradient_must_die_for_neg2"]:
        record(abs(sg) < 0.3,
               "THE MECHANISM: the baseline gradient is a property of composition-matched "
               "negatives and DIES for other-RBPs'-sites negatives", f"{sg:+.3f}", "|rho| < 0.3")

    rk = must("datasets where the dinuc baseline is lower than the GC baseline")
    if rk is not None:
        record(int(rk) == spec["n_rank_confounded"],
               "dn vs gc is rank-confounded with the baseline in every dataset", int(rk),
               spec["n_rank_confounded"])


def verify_baseline_order(T, g):
    """R1o: most of the magnitude is where the baseline stops; the fold range is not."""
    print("\nR1o baseline order  (is the 'contribution' just one order of composition up?)")
    d = T.get("baseline_order.csv")
    if d is None:
        return record(False, "baseline_order.csv present", "MISSING",
                      "run scripts/baseline_order.py")
    spec = g["r1o_baseline_order"]
    q = d.set_index("check")

    def must(k):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, "value"])

    n = q["n"].iloc[0] if "n" in q.columns else None
    if n is not None:
        record(int(n) == spec["n_datasets"], "datasets in the subsample", int(n),
               spec["n_datasets"])
    # THE ANCHOR, FIRST. Order 2 IS the paper's composition baseline, so it must reproduce the
    # paper's gain cell by cell. The first version of this script re-implemented the baseline,
    # overstated neg2's gain 2.31x, and agreed with the published R1 CONTRAST to 1e-4 by two
    # offsetting errors. An aggregate agreeing is not the baseline being the same baseline.
    rep = must("max |order-2 gain - published gain|")
    if rep is not None:
        at_most("order 2 reproduces the published per-dataset gain",
                rep, spec["max_order2_reproduction"])
    # KEYS ARE SPELLED OUT LITERALLY. test_golden_keys_are_read.py reads this file as TEXT, so
    # spec[f"gain_order{order}_{arm}"] is invisible to it and the key reads as unread. That is
    # tests/unit/test_golden_keys_are_read.py reads this file as text, so a key built with an
    # f-string reads as unread. Spell keys out literally.
    for label, key in (
            ("gain over order-2 baseline, gc arm", spec["gain_order2_gc"]),
            ("gain over order-3 baseline, gc arm", spec["gain_order3_gc"]),
            ("gain over order-2 baseline, dn arm", spec["gain_order2_dn"]),
            ("gain over order-3 baseline, dn arm", spec["gain_order3_dn"]),
            ("gain over order-2 baseline, neg2 arm", spec["gain_order2_neg2"]),
            ("gain over order-3 baseline, neg2 arm", spec["gain_order3_neg2"]),
            ("fraction removed by order 3, gc arm", spec["removed_fraction_gc"]),
            ("fraction removed by order 3, dn arm", spec["removed_fraction_dn"])):
        v = must(label)
        if v is not None:
            near(label, v, key)

    c2 = must("R1 contrast (dn-gc), order-2 baseline")
    c3 = must("R1 contrast (dn-gc), order-3 baseline")
    if c2 is not None:
        near("R1 contrast at order 2 (the subsample's own)", c2, spec["contrast_order2"])
        # the subsample must reproduce the full panel, or nothing here transfers
        pub = T.get("scale_check.csv")
        if pub is not None and "check" in pub.columns:
            pq = pub.set_index("check")
            k = "CONTRAST, AUROC scale (published headline)"
            if k in pq.index:
                at_most("the subsample reproduces the published R1 contrast",
                        abs(c2 - float(pq.loc[k, "value"])),
                        spec["max_subsample_contrast_drift"])
    if c3 is not None:
        near("R1 contrast at order 3", c3, spec["contrast_order3"])
    rf = must("fraction of the R1 contrast removed by order 3")
    if rf is not None:
        near("fraction of the R1 contrast removed by order 3", rf,
             spec["removed_fraction_contrast"])
        if spec["most_of_the_magnitude_must_be_removable"]:
            record(rf > 0.4, "MOST of the contrast's magnitude is removable by extending the "
                             "baseline one order, so 'sequence model' must be qualified",
                   f"{rf:.0%}", "> 40%")

    f2 = must("fold range across protocols, order-2 baseline")
    f3 = must("fold range across protocols, order-3 baseline")
    if f2 is not None:
        near("fold range at order 2", f2, spec["fold_range_order2"])
    if f3 is not None:
        near("fold range at order 3", f3, spec["fold_range_order3"])
    if f3 is not None:
        # The claim is "never collapses", NOT "widens". At n = 30 the change in the range is
        # [-1.89, +16.43] with P(<=0) = 0.18, so its direction is not resolved; an earlier
        # version of this gate asserted the range was unchanged, and the version after that
        # asserted it widened. Both overstated what 30 datasets can say.
        record(f3 >= spec["min_fold_range_order3"],
               "the fold range does not collapse when the baseline rises one order, so the "
               "paper's claim is not an artefact of where the baseline stops",
               f"{f3:.2f}x", f">= {spec['min_fold_range_order3']}x")
    if spec["range_change_direction_unresolved"]:
        k = "change in fold range, order 3 minus order 2"
        if k in q.index:
            lo = q.loc[k, "ci_low"]
            record(not pd.isna(lo) and float(lo) < 0,
                   "the change in the fold range straddles zero, so the paper must NOT claim "
                   "the range widens with baseline order",
                   f"CI low {float(lo):+.2f}" if not pd.isna(lo) else "NaN", "< 0")
        else:
            record(False, "the change in the fold range is reported with an interval",
                   "MISSING", "the row")

    # THE PER-DATASET SIGN. A ratio of means cannot show that the typical neg2 dataset retains
    # nothing; the median per-dataset ratio there is NEGATIVE.
    for label, key in (
            ("order-3 gain positive in, gc arm", spec["order3_positive_gc"]),
            ("order-3 gain positive in, dn arm", spec["order3_positive_dn"]),
            ("order-3 gain positive in, neg2 arm", spec["order3_positive_neg2"]),
            ("median per-dataset order-3/order-2 ratio, neg2 arm", spec["median_ratio_neg2"]),
            ("spearman(dataset size, order-3 retained fraction), dn arm",
             spec["size_dependence_dn"])):
        v = must(label)
        if v is not None:
            near(label, v, key)
    if spec["neg2_median_must_be_reported"]:
        v = must("median per-dataset order-3/order-2 ratio, neg2 arm")
        if v is not None:
            record(v <= 0, "the TYPICAL neg2 dataset retains nothing over a trinucleotide "
                           "baseline, so the positive panel mean must not be read as "
                           "'small but positive everywhere'", f"median {v:+.3f}", "<= 0")
    if spec["size_dependence_must_be_reported"]:
        record("spearman(dataset size, order-3 retained fraction), gc arm" in q.index,
               "the order-3 retained fraction's size dependence is reported, so the removal "
               "fractions cannot be read as unbiased for the full panel",
               "size dependence", "present")


def verify_baseline_order_models(T, g):
    """R1r: an order-3 baseline absorbs a near-constant amount from every model class.

    Not "the collapse is the 4-mer's alone": the absolute absorption is near-identical across
    model classes, so the differing shares are a denominator effect.
    """
    print("\nR1r baseline order across model classes  (is the collapse general?)")
    d = T.get("baseline_order_models.csv")
    if d is None:
        return record(False, "baseline_order_models.csv present", "MISSING",
                      "run scripts/baseline_order_models.py")
    spec = g["r1r_baseline_order_models"]
    q = d.set_index("check")

    def must(k):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, "value"])

    n = q["n"].iloc[0] if "n" in q.columns else None
    if n is not None:
        record(int(n) == spec["n_datasets"], "datasets", int(n), spec["n_datasets"])
    rep = must("max |order-2 gain - R1g published gain|")
    if rep is not None:
        at_most("order 2 reproduces R1g's per-dataset gain", rep,
                spec["max_order2_reproduction"])

    surv = {}
    for label, key in (
            ("kmer fraction surviving order 3, gc arm", spec["surviving_kmer_gc"]),
            ("cnn fraction surviving order 3, gc arm", spec["surviving_cnn_gc"]),
            ("splicebert fraction surviving order 3, gc arm", spec["surviving_splicebert_gc"]),
            ("kmer fraction surviving order 3, dn arm", spec["surviving_kmer_dn"]),
            ("cnn fraction surviving order 3, dn arm", spec["surviving_cnn_dn"]),
            ("splicebert fraction surviving order 3, dn arm", spec["surviving_splicebert_dn"])):
        v = must(label)
        if v is not None:
            near(label, v, key)
            surv[label] = v

    # THE SHARE GAP, and read the label carefully: it pins the gap in the SHARES, which is a
    # denominator effect. It is NOT evidence of differential fragility -- the absolute
    # absorption is near-identical and is gated separately below.
    for arm in ("gc", "dn"):
        k = surv.get(f"kmer fraction surviving order 3, {arm} arm")
        s = surv.get(f"splicebert fraction surviving order 3, {arm} arm")
        if k is not None and s is not None:
            record(s - k >= spec["min_deep_advantage_over_kmer"],
                   f"SpliceBERT's SHARE surviving the order-3 baseline exceeds the 4-mer's, "
                   f"{arm} arm (a denominator effect, not fragility -- see the absolute "
                   f"absorption)", f"{s - k:+.3f}", f">= {spec['min_deep_advantage_over_kmer']}")

    for label, key in (
            ("kmer gain over order-3 baseline, gc arm", spec["gain_order3_kmer_gc"]),
            ("splicebert gain over order-3 baseline, gc arm", spec["gain_order3_splicebert_gc"]),
            ("kmer gain over order-3 baseline, dn arm", spec["gain_order3_kmer_dn"]),
            ("splicebert gain over order-3 baseline, dn arm", spec["gain_order3_splicebert_dn"]),
            ("kmer R1 contrast (dn-gc), order-3 baseline", spec["contrast_order3_kmer"]),
            ("cnn R1 contrast (dn-gc), order-3 baseline", spec["contrast_order3_cnn"]),
            ("splicebert R1 contrast (dn-gc), order-3 baseline",
             spec["contrast_order3_splicebert"]),
            ("kmer dn/gc multiplier, order-3 baseline", spec["multiplier_order3_kmer"]),
            ("cnn dn/gc multiplier, order-3 baseline", spec["multiplier_order3_cnn"]),
            ("splicebert dn/gc multiplier, order-3 baseline",
             spec["multiplier_order3_splicebert"])):
        v = must(label)
        if v is not None:
            near(label, v, key)

    # THE PROTOCOL CONTRAST MUST SURVIVE AT ORDER 3 FOR EVERY MODEL, two-sided. A one-sided
    # floor: ">= x" would pass a rewritten value that is far too large.
    for model in ("kmer", "cnn", "splicebert"):
        lo = q.loc[f"{model} R1 contrast (dn-gc), order-3 baseline", "ci_low"] \
            if f"{model} R1 contrast (dn-gc), order-3 baseline" in q.index else None
        if lo is not None and not pd.isna(lo):
            record(float(lo) > 0,
                   f"protocol dependence survives an order-3 baseline for the {model}",
                   f"CI low {float(lo):+.4f}", "> 0")

    # THE ABSOLUTE COLUMN, which is what makes the shares readable. R1r's first version led
    # with the shares and drew a fragility conclusion the absolute numbers refute.
    absorbed = {}
    for label, key in (
            ("kmer absolute absorbed by order 3, gc arm", spec["absorbed_kmer_gc"]),
            ("cnn absolute absorbed by order 3, gc arm", spec["absorbed_cnn_gc"]),
            ("splicebert absolute absorbed by order 3, gc arm", spec["absorbed_splicebert_gc"]),
            ("kmer absolute absorbed by order 3, dn arm", spec["absorbed_kmer_dn"]),
            ("splicebert absolute absorbed by order 3, dn arm", spec["absorbed_splicebert_dn"])):
        v = must(label)
        if v is not None:
            near(label, v, key)
            absorbed[label] = v
    gcv = [v for k, v in absorbed.items() if k.endswith("gc arm")]
    if len(gcv) == 3:
        at_most("the order-3 baseline absorbs a near-CONSTANT absolute amount across model "
                "classes, so the differing shares are a denominator effect and not differing "
                "fragility", max(gcv) - min(gcv), spec["max_absorbed_spread_gc"])
    # The dinucleotide arm's spread. Gating the GC arm alone left this unchecked, and it is
    # roughly twice as large.
    dnv = must("absorbed spread across model classes, dn arm")
    if dnv is not None:
        at_most("and the same holds on the dinucleotide arm", dnv,
                spec["max_absorbed_spread_dn"])

    # THE COMPRESSION CORRECTION. The order-3 baseline raises composition, so part of the
    # "absorbed" amount is the moving ceiling -- and it takes MORE from a model with a bigger
    # gain, which is what made the raw absorption look constant. Corrected, the k-mer really
    # does lose more, by about 1.3-1.4x rather than the 3.4x the shares imply. BOTH of this
    # section's earlier framings were partial and the corrected ratio is the honest one.
    for label, key in (
            ("kmer compression-corrected order-3 residual, gc arm", spec["residual_kmer_gc"]),
            ("cnn compression-corrected order-3 residual, gc arm", spec["residual_cnn_gc"]),
            ("splicebert compression-corrected order-3 residual, gc arm",
             spec["residual_splicebert_gc"]),
            ("kmer compression-corrected order-3 residual, dn arm", spec["residual_kmer_dn"]),
            ("splicebert compression-corrected order-3 residual, dn arm",
             spec["residual_splicebert_dn"]),
            ("kmer/splicebert corrected-residual ratio, gc arm", spec["corrected_ratio_gc"]),
            ("kmer/splicebert corrected-residual ratio, dn arm", spec["corrected_ratio_dn"]),
            ("cnn minus kmer, absolute absorbed, gc arm", spec["cnn_minus_kmer_absorbed_gc"]),
            ("cnn minus kmer, absolute absorbed, dn arm", spec["cnn_minus_kmer_absorbed_dn"])):
        v = must(label)
        if v is not None:
            near(label, v, key)
    if spec["compression_correction_must_be_reported"]:
        record("kmer/splicebert corrected-residual ratio, gc arm" in q.index,
               "the compression-corrected comparison is reported, so neither the shares nor "
               "the raw absorption can be quoted as the whole story",
               "corrected residual ratio", "present")
        # AND THE CNN MUST STILL ABSORB THE LEAST, on the corrected scale too. Any capacity
        # story predicts otherwise, which is why this is the gate that kills capacity readings.
        rk = must("kmer compression-corrected order-3 residual, gc arm")
        rc = must("cnn compression-corrected order-3 residual, gc arm")
        rs = must("splicebert compression-corrected order-3 residual, gc arm")
        if None not in (rk, rc, rs):
            record(abs(rc) < min(abs(rk), abs(rs)),
                   "the CNN loses the LEAST even after the compression correction, which no "
                   "capacity story predicts",
                   f"cnn {rc:+.4f} vs kmer {rk:+.4f}, splicebert {rs:+.4f}", "cnn smallest")

    # THE PER-DATASET SIGN. A mean of +0.0058 is consistent with "small everywhere" and with
    # "negative in a third of the panel". For the 4-mer over a trinucleotide baseline it is
    # the second, and that is the result the shares could not show.
    pos = {}
    for label, key in (
            ("kmer order-3 gain positive in, gc arm", spec["kmer_order3_positive_gc"]),
            ("cnn order-3 gain positive in, gc arm", spec["cnn_order3_positive_gc"]),
            ("splicebert order-3 gain positive in, gc arm",
             spec["splicebert_order3_positive_gc"]),
            ("splicebert order-3 gain positive in, dn arm",
             spec["splicebert_order3_positive_dn"])):
        v = must(label)
        if v is not None:
            near(label, v, key)
            pos[label] = v
    if spec["splicebert_must_be_positive_everywhere"]:
        for arm in ("gc", "dn"):
            v = pos.get(f"splicebert order-3 gain positive in, {arm} arm")
            k = pos.get(f"kmer order-3 gain positive in, {arm} arm")
            if v is not None:
                record(int(v) == spec["n_datasets"],
                       f"SpliceBERT's contribution over a TRINUCLEOTIDE baseline is positive "
                       f"on every dataset, {arm} arm", f"{int(v)}/{spec['n_datasets']}",
                       f"{spec['n_datasets']}/{spec['n_datasets']}")
            if v is not None and k is not None:
                record(v > k, f"and the 4-mer's is not, {arm} arm",
                       f"kmer {int(k)} vs splicebert {int(v)}", "kmer fewer")

    if spec["splicebert_multiplier_must_be_lowest"]:
        m = {k: must(f"{k} dn/gc multiplier, order-3 baseline")
             for k in ("kmer", "cnn", "splicebert")}
        if all(v is not None for v in m.values()):
            record(m["splicebert"] < min(m["kmer"], m["cnn"]),
                   "SpliceBERT still has the LOWEST multiplier at order 3, so 'the contrast "
                   "grows with capacity' stays withdrawn",
                   f"{m['splicebert']:.2f}x", f"< {min(m['kmer'], m['cnn']):.2f}x")


def verify_models_by_protocol(T, g):
    """R1w: the floor and the recommendation generalise across models. R1n does not."""
    print("\nR1w the 4-mer-only analyses, per model class  (what generalises?)")
    d = T.get("models_by_protocol.csv")
    if d is None:
        return record(False, "models_by_protocol.csv present", "MISSING",
                      "run scripts/models_by_protocol.py")
    spec = g["r1w_models_by_protocol"]
    q = d.set_index("check")

    def must(k):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, "value"])

    for label, key in (
            ("minimum fold range over transforms, kmer", spec["floor_kmer"]),
            ("minimum fold range over transforms, cnn", spec["floor_cnn"]),
            ("minimum fold range over transforms, splicebert", spec["floor_splicebert"]),
            ("equalising exponent, kmer", spec["exponent_kmer"]),
            ("equalising exponent, cnn", spec["exponent_cnn"]),
            ("equalising exponent, splicebert", spec["exponent_splicebert"]),
            ("fold range at the equalising exponent, cnn", spec["argmin_range_cnn"]),
            ("fold range at the equalising exponent, splicebert",
             spec["argmin_range_splicebert"]),
            ("incremental R2 of the protocol label, kmer", spec["inc_protocol_kmer"]),
            ("incremental R2 of the baseline, kmer", spec["inc_baseline_kmer"]),
            ("incremental R2 of the protocol label, cnn", spec["inc_protocol_cnn"]),
            ("incremental R2 of the baseline, cnn", spec["inc_baseline_cnn"]),
            ("incremental R2 of the protocol label, splicebert",
             spec["inc_protocol_splicebert"]),
            ("incremental R2 of the baseline, splicebert", spec["inc_baseline_splicebert"]),
            ("within-arm spearman(baseline, gain), neg2, splicebert",
             spec["gradient_neg2_splicebert"]),
            ("within-arm spearman(baseline, gain), neg2, kmer", spec["gradient_neg2_kmer"])):
        v = must(label)
        if v is not None:
            near(label, v, key)

    # 1. No model class reaches protocol independence among the fixed coordinates. Gated on
    # the WEAKEST, because "no rescaling works" is a claim about all of them.
    floors = [must(f"minimum fold range over transforms, {m}")
              for m in ("kmer", "cnn", "splicebert")]
    if all(f is not None for f in floors):
        record(min(floors) >= spec["min_floor_over_models"],
               "every model class stays above the equal-means null among fixed coordinates, "
               "so the central negative result is not a 4-mer artefact",
               f"{min(floors):.2f}x (weakest)", f">= {spec['min_floor_over_models']}x")

    # 2. THE CNN REVERSAL MUST BE REPORTED. For the CNN the protocol LABEL beats the baseline,
    # so R1n's headline is false for it and the pooled 4-mer figure would conceal that.
    if spec["cnn_reversal_must_be_reported"]:
        pc, bc = (must("incremental R2 of the protocol label, cnn"),
                  must("incremental R2 of the baseline, cnn"))
        if pc is not None and bc is not None:
            record(pc > bc, "the CNN reversal is present and reported: for that model the "
                            "protocol LABEL explains more than the baseline, so 'it is the "
                            "baseline, not the label' is 4-mer-specific",
                   f"protocol {pc:.4f} > baseline {bc:.4f}", "reversed")

    # 3. AND THE TWO-FAMILY MECHANISM IS 4-MER-SPECIFIC. SpliceBERT's bias-aware gradient is
    # strong, where the k-mer's is absent; the mechanism cannot be stated as general.
    if spec["two_family_mechanism_is_kmer_specific"]:
        gk, gs = (must("within-arm spearman(baseline, gain), neg2, kmer"),
                  must("within-arm spearman(baseline, gain), neg2, splicebert"))
        if gk is not None and gs is not None:
            record(abs(gs) > abs(gk) + 0.2,
                   "SpliceBERT's bias-aware gradient is strong where the 4-mer's is absent, so "
                   "the two-family mechanism must be labelled a 4-mer result",
                   f"splicebert {gs:+.3f} vs kmer {gk:+.3f}", "much stronger")

    # 4. The recommendation does generalise, in all three protocol pairs for every model.
    for m in ("kmer", "cnn", "splicebert"):
        v = must(f"protocol pairs where rank agreement improves, {m}")
        if v is not None:
            record(int(v) == spec["recommendation_pairs_improved"],
                   f"the recommendation improves rank agreement in every protocol pair, {m}",
                   f"{int(v)}/3", "3/3")


def verify_three_arm_models(T, g):
    """R1v: the three-protocol span holds for every model class, and is widest for the CNN."""
    print("\nR1v three protocols x three models  (is the span a k-mer artefact?)")
    d = T.get("three_arm_models.csv")
    if d is None:
        return record(False, "three_arm_models.csv present", "MISSING",
                      "run scripts/three_arm_models.py")
    spec = g["r1v_three_arm_models"]
    q = d.set_index("check")

    def must(k):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, "value"])

    n = q["n"].iloc[0] if "n" in q.columns else None
    if n is not None:
        record(int(n) == spec["n_datasets"], "datasets", int(n), spec["n_datasets"])
    for label, key in (
            ("protocol range of panel means, within kmer", spec["range_protocol_kmer"]),
            ("protocol range of panel means, within cnn", spec["range_protocol_cnn"]),
            ("protocol range of panel means, within splicebert",
             spec["range_protocol_splicebert"]),
            ("model-class range of panel means, within gc arm", spec["range_modelclass_gc"]),
            ("model-class range of panel means, within dn arm", spec["range_modelclass_dn"]),
            ("model-class range of panel means, within neg2 arm",
             spec["range_modelclass_neg2"]),
            ("composition AUROC, dn arm", spec["comp_dn"]),
            ("composition AUROC, gc arm", spec["comp_gc"]),
            ("composition AUROC, neg2 arm", spec["comp_neg2"]),
            ("kmer nested contribution, neg2 arm", spec["kmer_neg2"]),
            ("cnn nested contribution, neg2 arm", spec["cnn_neg2"]),
            ("splicebert nested contribution, neg2 arm", spec["splicebert_neg2"]),
            ("kmer three-protocol span", spec["span_kmer"]),
            ("cnn three-protocol span", spec["span_cnn"]),
            ("splicebert three-protocol span", spec["span_splicebert"])):
        v = must(label)
        if v is not None:
            near(label, v, key)

    # THE CLAIM. Gated on the WEAKEST-spanning model, because the sentence the paper makes is
    # "for every model class", and a mean over models would let one collapse unnoticed.
    lo = must("minimum three-protocol span over model classes")
    if lo is not None:
        record(lo >= spec["min_span_over_models"],
               "the three-protocol span survives for EVERY model class, so it is not an "
               "artefact of measuring with a bag of k-mers",
               f"{lo:.2f}x (weakest)", f">= {spec['min_span_over_models']}x")

    # THE REVERSAL, generalised. The bias-aware protocol is the easiest discrimination and
    # yields the least for every model, which is the paper's title claim.
    if spec["neg2_lowest_for_all_models"]:
        v = must("model classes whose lowest contribution is the bias-aware arm")
        if v is not None:
            record(int(v) == 3, "the bias-aware arm yields the LEAST for every model class, "
                                "which is the reversal the title states", f"{int(v)}/3", "3/3")

    # AND IT DOES NOT GROW WITH CAPACITY. The CNN spans widest, not SpliceBERT. If that ever
    # inverts, "the contrast grows with capacity" would need re-examining rather than the
    # withdrawal being quietly retained.
    if spec["span_must_not_grow_with_capacity"]:
        k, c, s = (must("kmer three-protocol span"), must("cnn three-protocol span"),
                   must("splicebert three-protocol span"))
        if None not in (k, c, s):
            record(s < c and s < k,
                   "SpliceBERT has the NARROWEST span of the three, so the span does not grow "
                   "with capacity and the withdrawn claim stays withdrawn",
                   f"splicebert {s:.2f}x vs cnn {c:.2f}x, kmer {k:.2f}x", "splicebert lowest")

    # A5's SURVIVING CLAIM. In absolute AUROC neither effect dominates -- which is larger
    # depends on the arm and the model -- so the only ordering the paper may assert is the one
    # in fold terms: the protocol span within a model exceeds the model-class span within an
    # arm. Read off the note column, where each range carries its own fold, so a reader cannot
    # mix an AUROC difference with a fold change the way the Introduction once did.
    folds = {}
    for r in q.index:
        if str(r).startswith(("protocol range of panel means",
                              "model-class range of panel means")):
            note = str(q.loc[r, "note"])
            if note.startswith("fold "):
                folds[str(r)] = float(note.split()[1])
    prot = [v for k_, v in folds.items() if k_.startswith("protocol")]
    mcls = [v for k_, v in folds.items() if k_.startswith("model-class")]
    if len(prot) == 3 and len(mcls) == 3 and spec["fold_ranges_must_overlap"]:
        # THE FIRST VERSION OF THIS GATE ASSERTED THAT THE PROTOCOL SPAN DOMINATES, AND FAILED.
        # It should have: protocol folds run 3.76 to 7.63 and model-class folds 2.65 to 4.14,
        # so SpliceBERT's protocol span sits BELOW the bias-aware arm's model-class span and
        # the two ranges overlap. What is true is that the protocol range reaches higher while
        # still overlapping, i.e. the effects are comparable and neither dominates. That is
        # what the paper may say, so it is what is gated.
        record(max(prot) > max(mcls) and min(prot) < max(mcls),
               "the protocol and model-class fold ranges OVERLAP while the protocol range "
               "reaches higher, so the two effects are comparable and neither dominates",
               f"protocol {min(prot):.2f}-{max(prot):.2f}x, "
               f"model class {min(mcls):.2f}-{max(mcls):.2f}x",
               "overlapping, protocol reaching higher")


def verify_transport(T, g):
    """R1u: an equalising exponent exists in sample and does not transport. Both halves gated."""
    print("\nR1u transportability  (does any rescaling carry to another benchmark?)")
    d = T.get("transport_check.csv")
    if d is None:
        return record(False, "transport_check.csv present", "MISSING",
                      "run scripts/transport_check.py")
    spec = g["r1u_transport"]
    q = d.set_index("check")

    def must(k):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, "value"])

    k = "external rank agreement, raw"
    if k in q.index and "n" in d.columns:
        record(int(q.loc[k, "n"]) == spec["n_datasets_external"],
               "the external tests run on Horlacher's full shared panel",
               int(q.loc[k, "n"]), spec["n_datasets_external"])
    for label, key in (
            ("fold range at the headroom coordinate, p = 1", spec["range_at_p1"]),
            ("equalising exponent on our protocols", spec["exponent_ours"]),
            ("fold range at our equalising exponent", spec["range_at_argmin"]),
            ("equalising exponent on Horlacher's protocols", spec["exponent_theirs"]),
            ("their fold range, raw", spec["their_range_raw"]),
            ("their fold range at OUR exponent", spec["their_range_at_our_p"]),
            ("their fold range at their own exponent", spec["their_range_at_their_p"]),
            ("equal-means null fold range, median", spec["null_median"]),
            ("equal-means null fold range, 95th percentile", spec["null_p95"]),
            ("external rank agreement, raw", spec["external_rank_raw"]),
            ("external rank agreement, headroom", spec["external_rank_headroom"]),
            ("external scale-free disagreement, raw", spec["external_disagreement_raw"]),
            ("external scale-free disagreement, headroom",
             spec["external_disagreement_headroom"]),
            ("external within-arm spearman(baseline, gain), n1",
             spec["external_gradient_n1"]),
            ("external within-arm spearman(baseline, gain), n2",
             spec["external_gradient_n2"])):
        v = must(label)
        if v is not None:
            near(label, v, key)

    # 1. THE IN-SAMPLE EQUALISING EXPONENT MUST BE REPORTED. Suppressing it restores a claim
    # that is false in the paper's own data, which is the worst kind of claim to hold.
    if spec["in_sample_exponent_must_be_reported"]:
        v = must("fold range at our equalising exponent")
        if v is not None:
            record(v < 1.1, "an in-sample equalising exponent EXISTS and is reported, so the "
                            "paper cannot claim that no rescaling reaches a protocol-free "
                            "quantity", f"{v:.3f}x at p = 1.544", "< 1.1x")

    # 1b. AND IT MUST BE SHOWN TO BE AN AGGREGATION ARTEFACT. g/(1-c)^p is a ratio and 1-c
    # reaches 0.0267 here, so a MEAN of that ratio is dominated by a few high-baseline cells
    # and increasingly so with p. That mean is the published aggregation and it is the only one
    # that touches 1. Both halves are gated: the mean-of-ratios minimum must still fall below
    # 1.1, so the inconvenient in-sample fact cannot be suppressed, AND the two aggregations
    # that are not a mean of a ratio must stay clear of 1, which is what makes it an artefact
    # rather than a coordinate. If they ever converge, the retraction has to be revisited.
    mean_min = must("minimum span over p, mean of per-dataset ratios")
    med_min = must("minimum span over p, median of per-dataset ratios")
    rom_min = must("minimum span over p, ratio of panel means")
    for label, key in (
            ("minimum span over p, median of per-dataset ratios", spec["min_span_median"]),
            ("minimum span over p, ratio of panel means", spec["min_span_ratio_of_means"]),
            ("min headroom 1-c over all cells", spec["min_headroom"])):
        v = must(label)
        if v is not None:
            near(label, v, key)
    n_low = must("cells with headroom below 0.15")
    if n_low is not None:
        at_least("cells whose headroom is small enough to dominate a mean of a ratio",
                 int(n_low), spec["min_cells_below_headroom_015"])
    if None not in (mean_min, med_min, rom_min) and spec["aggregation_artefact_must_hold"]:
        floor = spec["min_span_under_robust_aggregation"]
        record(mean_min < 1.1 and med_min > floor and rom_min > floor,
               "the equalising exponent closes the span ONLY under a mean of ratios, so it is "
               "an aggregation artefact and the paper's retraction stands",
               f"mean {mean_min:.3f}, median {med_min:.3f}, ratio-of-means {rom_min:.3f}",
               f"mean < 1.1 and both others > {floor}")

    # 2. AND IT MUST NOT TRANSPORT. This is now the paper's claim, so it is the thing that has
    # to be checked -- not a floor.
    pa, pt = must("equalising exponent on our protocols"), \
        must("equalising exponent on Horlacher's protocols")
    if pa and pt:
        record(pt / pa >= spec["min_exponent_ratio"],
               "the equalising exponent does NOT transport between benchmarks, which is the "
               "claim that replaces the retracted floor", f"{pt / pa:.2f}x apart",
               f">= {spec['min_exponent_ratio']}x")
    ours_there = must("their fold range at OUR exponent")
    raw_there = must("their fold range, raw")
    if ours_there and raw_there:
        record(ours_there > 0.9 * raw_there,
               "and our exponent buys essentially nothing on their benchmark, so the "
               "normalisation is fitted rather than universal",
               f"{ours_there:.2f}x vs raw {raw_there:.2f}x", "no material gain")

    # 3. THE 2.00x INFERENCE, against the CORRECT null. max/min over three noisy means is
    # bounded below by 1, so comparing against 1.0 was too generous.
    obs, p95 = must("fold range at the headroom coordinate, p = 1"), \
        must("equal-means null fold range, 95th percentile")
    if obs and p95:
        record(obs > p95, "the headroom coordinate's range exceeds the EQUAL-MEANS null's 95th "
                          "percentile, not merely 1.0, so the range is not an artefact of "
                          "taking a ratio of three noisy means",
               f"{obs:.3f}x vs null p95 {p95:.3f}x", "outside the null")

    # 4. THE EXTERNAL FALSIFICATION MUST BE REPORTED whichever way it points. It points against
    # the recommendation, and the paper says so.
    if spec["external_falsification_must_be_reported"]:
        for k in ("external rank agreement, headroom",
                  "external scale-free disagreement, headroom"):
            record(k in q.index, f"the external test of the recommendation is reported ({k})",
                   k, "present")
        rr, rh = must("external rank agreement, raw"), must("external rank agreement, headroom")
        if rr is not None and rh is not None:
            record(rh < rr, "and it FAILS to replicate out of sample -- rank agreement falls "
                            "under normalisation, which is this script's own pre-registered "
                            "falsification criterion", f"{rr:+.3f} -> {rh:+.3f}", "falls")


def verify_fold_integrity(T, g):
    """Did the committed scores use the study's chromosome folds? For 20 of 188, no."""
    print("\nfold integrity  (did the committed scores use config/folds.tsv?)")
    d = T.get("fold_integrity.csv")
    if d is None:
        return record(False, "fold_integrity.csv present", "MISSING",
                      "run scripts/fold_integrity.py")
    spec = g["fold_integrity"]
    q = d.set_index("check")

    def must(k):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, "value"])

    total = ((must("score sets audited, gc arm") or 0)
             + (must("score sets audited, dn arm") or 0)
             + (must("score sets audited, neg2 arm") or 0))
    record(int(total) == spec["n_score_sets"], "(dataset, arm, model) score sets audited",
           int(total), spec["n_score_sets"])

    # Keys spelled out literally -- test_golden_keys_are_read.py reads this file as TEXT and
    # cannot see an f-string, so a key built by interpolation reads as unread.
    for label, key in (
            ("datasets NOT chromosome-grouped, gc arm",
             {"value": 94 - spec["grouped_gc"]["value"], "tol": 0}),
            ("datasets NOT chromosome-grouped, dn arm",
             {"value": 94 - spec["grouped_dn"]["value"], "tol": 0}),
            ("datasets NOT chromosome-grouped, neg2 arm",
             {"value": 94 - spec["grouped_neg2"]["value"], "tol": 0}),
            ("datasets aligned to folds.tsv, gc arm", spec["aligned_gc"]),
            ("datasets aligned to folds.tsv, neg2 arm", spec["aligned_neg2"]),
            ("max chromosomes in any score fold, gc arm",
             spec["max_chroms_per_score_fold_gc"]),
            ("max chromosomes in any score fold, dn arm",
             spec["max_chroms_per_score_fold_dn"]),
            ("max chromosomes in any score fold, neg2 arm",
             spec["max_chroms_per_score_fold_neg2"]),
            ("leaky datasets", spec["leaky_datasets"]),
            ("datasets aligned to folds.tsv, dn arm", spec["aligned_dn"])):
        v = must(label)
        if v is not None:
            near(label, v, key)

    # 1. THE GC ARM CARRIES THE HEADLINE. It must be perfectly chromosome-grouped, and the
    # direct leakage metric must be exactly zero -- not small, zero, because under chromosome
    # grouping a within-1kb neighbour is necessarily in the same fold by construction.
    if spec["gc_must_be_fully_grouped"]:
        n_bad = must("datasets NOT chromosome-grouped, gc arm")
        if n_bad is not None:
            record(int(n_bad) == 0,
                   "the GC arm, which carries every headline, is fully chromosome-grouped",
                   f"{int(n_bad)} leaky", "0")
        nbr = must("max cross-fold 1kb neighbour fraction, gc arm")
        if nbr is not None:
            at_most("and its direct cross-fold neighbour fraction is exactly zero, as "
                    "chromosome grouping requires by construction",
                    nbr, spec["max_cross_fold_neighbours_gc"])

    # 1b. AND THE BIAS-AWARE ARM, to the same standard. It is the denominator of every span,
    # and until this gate existed the manuscript's claim that it was clean rested on an
    # argument from construction rather than on a measurement.
    if spec["neg2_must_be_fully_grouped"]:
        n_bad2 = must("datasets NOT chromosome-grouped, neg2 arm")
        if n_bad2 is not None:
            record(int(n_bad2) == 0,
                   "the bias-aware arm, the denominator of every span, is fully "
                   "chromosome-grouped", f"{int(n_bad2)} leaky", "0")
        nbr2 = must("max cross-fold 1kb neighbour fraction, neg2 arm")
        if nbr2 is not None:
            at_most("and its direct cross-fold neighbour fraction is exactly zero too",
                    nbr2, spec["max_cross_fold_neighbours_neg2"])

    # 1c. AND THE DINUCLEOTIDE ARM, which is where the defect was. It used to be held to a
    # weaker standard -- 74 of 94 grouped, 20 leaky, with a bound on how far the conclusion
    # could move when the 20 were dropped -- because the sweep had not been rerun. It has been,
    # so the arm is held to the same zero as the other two and the tolerance is gone.
    if spec["dn_must_be_fully_grouped"]:
        n_bad3 = must("datasets NOT chromosome-grouped, dn arm")
        if n_bad3 is not None:
            record(int(n_bad3) == 0,
                   "the dinucleotide arm, where the partition defect was, is fully "
                   "chromosome-grouped after the retrain", f"{int(n_bad3)} leaky", "0")
        nbr3 = must("max cross-fold 1kb neighbour fraction, dn arm")
        if nbr3 is not None:
            at_most("and its direct cross-fold neighbour fraction is exactly zero, which it "
                    "was not before: it reached 44.5%",
                    nbr3, spec["max_cross_fold_neighbours_dn"])

    # 2. NO DATASET IN ANY ARM MAY BE LEAKY. Pinned two-sided at zero: a floor would pass a
    # rerun that leaked more, and the whole point of the retrain is that zero is now reachable.
    nl = must("leaky datasets")
    if nl is not None:
        record(int(nl) == spec["leaky_datasets"]["value"],
               "no dataset in any arm is off-partition, so the defect is corrected rather "
               "than disclosed", int(nl), spec["leaky_datasets"]["value"])


def verify_estimands(T, g):
    """B1: is the protocol ordering a property of AUROC, or of the comparison?"""
    print("\nestimands  (does the ordering survive leaving the ROC?)")
    d = T.get("estimands.csv")
    if d is None:
        return record(False, "estimands.csv present", "MISSING",
                      "run scripts/estimands.py --store ../rbp-store")
    spec = g["estimands"]
    q = d.set_index("check")

    def get(k):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, "value"])

    # THE HEADLINE. Five estimands, three model classes, one ordering.
    n = get("estimands on which the protocol ordering holds for every model")
    if n is not None:
        record(int(n) == spec["n_estimands_ordered"],
               "dinucleotide > GC > bias-aware on EVERY estimand for EVERY model class, so "
               "the direction of the protocol effect is not an artefact of measuring in AUROC",
               f"{int(n)}/5", f"{spec['n_estimands_ordered']}/5")

    # AND THE MAGNITUDE, which is scale-dependent and must not be quoted as if it were not.
    spans = {}
    for model, key in (("kmer", spec["deviance_span_kmer"]),
                       ("cnn", spec["deviance_span_cnn"]),
                       ("splicebert", spec["deviance_span_splicebert"])):
        v = [get(f"delta_deviance, {a} arm, {model}") for a in ("dn", "gc", "neg2")]
        if None in v or min(v) <= 0:
            continue
        spans[model] = max(v) / min(v)
        near(f"delta-deviance span, {model}", spans[model], key)
    if spans:
        at_least("the deviance span stays clear of 1, so 'roughly two-fold on an unbounded "
                 "scale' is a smaller effect than the AUROC figure and not an absent one",
                 min(spans.values()), spec["min_deviance_span"])

    # THE RESIDUAL DIAGNOSTIC. A small increment is not an absent signal.
    res = [float(q.loc[k, "value"]) for k in q.index
           if str(k).startswith("residual_auroc,")]
    if res:
        near("residual AUROC, kmer, dinucleotide arm",
             get("residual_auroc, dn arm, kmer") or float("nan"),
             spec["residual_auroc_kmer_dn"])
        at_least("orthogonalising every model score against the composition block still "
                 "leaves discrimination above chance in every arm, so a small increment does "
                 "not mean an absent signal", min(res),
                 spec["residual_auroc_must_exceed_chance"])


def verify_nested_scale(T, g):
    """1b/1c: the covariate scale and the standardisation window, both measured."""
    print("\nnested scale  (does the covariate scale or the standardisation window matter?)")
    d = T.get("nested_scale.csv")
    if d is None:
        return record(False, "nested_scale.csv present", "MISSING",
                      "run scripts/nested_scale.py --store ../rbp-store")
    spec = g["nested_scale"]
    q = d.set_index("check")

    def get(k):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, "value"])

    for label, key in (
            ("logit minus probability, gc arm, cnn", spec["logit_gain_cnn_gc"]),
            ("logit minus probability, gc arm, splicebert", spec["logit_gain_splicebert_gc"]),
            ("logit minus probability, dn arm, cnn", spec["logit_gain_cnn_dn"]),
            ("logit minus probability, dn arm, splicebert", spec["logit_gain_splicebert_dn"])):
        v = get(label)
        if v is not None:
            near(label, v, key)

    # THE DIRECTION IS THE CLAIM. "The published figures are conservative" is only true while
    # the logit scale raises every neural contribution in every arm.
    if spec["logit_must_raise_neural"]:
        rows = [k for k in q.index if str(k).startswith("logit minus probability")]
        vals = [float(q.loc[k, "value"]) for k in rows]
        record(len(vals) == 6 and all(v > 0 for v in vals),
               "the logit scale RAISES the neural contribution in all six arm-by-model cells, "
               "so the published probability-scale figures are the conservative ones",
               f"{sum(v > 0 for v in vals)}/{len(vals)} positive", "6/6")

    # THE CONTRAST IS WHAT THE PAPER IS ABOUT, so it gets its own ceiling. A scale change that
    # shifts both arms equally shifts no contrast, which is why the per-arm gains above are an
    # upper bound on this rather than a substitute for it.
    cs = [abs(float(q.loc[k, "value"])) for k in q.index
          if str(k).startswith("contrast shift under the logit scale")]
    if cs:
        at_most("the covariate scale moves the two-arm contrast in the fourth decimal, well "
                "inside the protein-clustered half-width, so the choice is not load-bearing",
                max(cs), spec["max_contrast_shift_from_scale"])

    # AND THE STANDARDISATION WINDOW MUST STAY BELOW THE QUOTED PRECISION.
    wf = [abs(float(q.loc[k, "value"])) for k in q.index
          if str(k).startswith("within-fold minus whole-dataset")]
    if wf:
        at_most("standardising within fold instead of over the dataset changes no panel mean "
                "beyond the precision the paper quotes, so the improper form costs nothing",
                max(wf), spec["max_within_fold_effect"])


def verify_design_effect(T, g):
    """The design effect the paper applies, against the one measured from the data."""
    print("\ndesign effect  (1.35 is applied; what does the data say?)")
    d = T.get("design_effect.csv")
    if d is None:
        return record(False, "design_effect.csv present", "MISSING",
                      "run scripts/design_effect.py")
    spec = g["design_effect"]
    q = d.set_index("check")

    def must(k):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, "value"])

    for label, key in (
            ("ratio_clustering_10kb", spec["clustering_10kb"]),
            ("ratio_clustering_100kb", spec["clustering_100kb"]),
            ("ratio_clustering_1000kb", spec["clustering_1000kb"]),
            ("ratio_fitting", spec["fitting"]),
            ("measured design effect, 10 kb blocks", spec["measured_product"]),
            ("datasets significant at design effect 1.00", spec["n_significant_unadjusted"]),
            ("datasets significant at design effect 1.35", spec["n_significant_applied"])):
        v = must(label)
        if v is not None:
            near(label, v, key)

    # THE CLAIM: the applied figure is the conservative one. If the measured product ever
    # rises above 1.35 the paper is understating its own uncertainty and the sentence in
    # Methods has to be rewritten, so this is gated rather than assumed.
    if spec["applied_must_exceed_measured"]:
        m = must("measured design effect, 10 kb blocks")
        if m is not None:
            record(m < 1.35, "the applied design effect is the conservative one",
                   f"measured {m:.3f}", "< 1.35")


def verify_standalone_auroc(T, g):
    """The model's own AUROC per arm, which Table 1 needs and no table held until now."""
    print("\nstandalone AUROC  (the model alone, not composition plus the model)")
    d = T.get("standalone_auroc.csv")
    if d is None:
        return record(False, "standalone_auroc.csv present", "MISSING",
                      "run scripts/standalone_auroc.py")
    spec = g["standalone_auroc"]
    q = d.set_index("check")

    def must(k):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, "value"])

    for label, key in (
            ("model alone, dn arm", spec["model_alone_dn"]),
            ("model alone, gc arm", spec["model_alone_gc"]),
            ("model alone, neg2 arm", spec["model_alone_neg2"]),
            ("model-alone AUROC drop, gc to dn", spec["drop_gc_to_dn"]),
            ("datasets where composition beats the model, dn arm", spec["comp_beats_model_dn"]),
            ("datasets where composition beats the model, gc arm", spec["comp_beats_model_gc"]),
            ("datasets where composition beats the model, neg2 arm",
             spec["comp_beats_model_neg2"])):
        v = must(label)
        if v is not None:
            near(label, v, key)

    for label, key in (
            ("pearson(model alone, composition alone), pooled over arms",
             spec["corr_alone_comp_pooled"]),
            ("pearson(model alone, composition alone), gc arm", spec["corr_alone_comp_gc"]),
            ("pearson(model alone, composition alone), dn arm", spec["corr_alone_comp_dn"]),
            ("pearson(model alone, composition alone), neg2 arm",
             spec["corr_alone_comp_neg2"])):
        v = must(label)
        if v is not None:
            near(label, v, key)
    k = "pearson(model alone, composition alone), pooled over arms"
    if k in q.index and spec["corr_alone_comp_must_exclude_zero"]:
        lo = float(q.loc[k, "ci_low"])
        record(lo > 0.5, "apparent difficulty and composition difficulty are the SAME axis, "
                         "which is what licenses treating the baseline as the channel the "
                         "protocol acts through", f"clustered lower bound {lo:.4f}", "> 0.5")

    # B7. THE 3x3. Values first, then the two orderings, which are what the text may say.
    for label, key in (
            ("model alone, kmer, dn arm", spec["alone_kmer_dn"]),
            ("model alone, cnn, dn arm", spec["alone_cnn_dn"]),
            ("model alone, splicebert, dn arm", spec["alone_splicebert_dn"]),
            ("model alone, kmer, neg2 arm", spec["alone_kmer_neg2"]),
            ("model alone, cnn, neg2 arm", spec["alone_cnn_neg2"]),
            ("model alone, splicebert, neg2 arm", spec["alone_splicebert_neg2"]),
            ("splicebert bias-aware minus GC, model alone",
             spec["splicebert_neg2_minus_gc"])):
        v = must(label)
        if v is not None:
            near(label, v, key)
    for label, key in (
            ("datasets where difficulty and contribution move oppositely, GC to dinucleotide",
             spec["opposite_gc_to_dn"]),
            ("datasets where difficulty and contribution move oppositely, GC to bias-aware",
             spec["opposite_gc_to_neg2"]),
            ("spearman(delta model-alone, delta gain), GC to dinucleotide",
             spec["magnitude_corr_gc_to_dn"])):
        v = must(label)
        if v is not None:
            near(label, v, key)

    nh, ne = (must("models for which dinucleotide matching is the hardest"),
              must("models for which the bias-aware arm is the easiest"))
    if nh is not None:
        record(int(nh) == spec["models_dn_hardest"],
               "dinucleotide matching is the hardest discrimination for EVERY model class, so "
               "that half of the inversion is not a property of one measurement",
               f"{int(nh)}/3", f"{spec['models_dn_hardest']}/3")
    if ne is not None:
        record(int(ne) == spec["models_neg2_easiest"],
               "the bias-aware arm is the easiest for only TWO of three model classes by their "
               "own AUROC, so 'the easiest protocol' needs the qualification the text gives it",
               f"{int(ne)}/3", f"{spec['models_neg2_easiest']}/3")

    # THE PROGRESSION IS THE CLAIM, and it is the sharpest form of the paper's point: the
    # easier a protocol looks, the more often nineteen composition features beat the model.
    if spec["comp_beats_model_must_increase_with_easiness"]:
        a = must("datasets where composition beats the model, dn arm")
        b = must("datasets where composition beats the model, gc arm")
        c = must("datasets where composition beats the model, neg2 arm")
        if None not in (a, b, c):
            record(a < b < c,
                   "composition beats the model more often as the protocol gets easier",
                   f"{int(a)} < {int(b)} < {int(c)}", "strictly increasing")


def verify_region_asymmetry(T, g):
    """Region is matched in two arms and free in the third. How much does that buy?"""
    print("\nregion asymmetry  (the bias-aware arm does not match transcript region)")
    d = T.get("region_asymmetry.csv")
    if d is None:
        return record(False, "region_asymmetry.csv present", "MISSING",
                      "run scripts/region_asymmetry.py")
    spec = g["region_asymmetry"]
    q = d.set_index("check")

    def must(k):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, "value"])

    for label, key in (
            ("region-only AUROC, gc arm", spec["region_auroc_gc"]),
            ("region-only AUROC, dn arm", spec["region_auroc_dn"]),
            ("region-only AUROC, neg2 arm", spec["region_auroc_neg2"]),
            ("region-only AUROC after matching, neg2 arm", spec["region_auroc_matched"]),
            ("composition alone, neg2 arm region-matched", spec["comp_matched"]),
            ("nested contribution, neg2 arm region-matched", spec["gain_matched"]),
            ("region-only AUROC, pipeline region-matched arm", spec["region_auroc_built"]),
            ("composition alone, pipeline region-matched arm", spec["comp_built"]),
            ("nested contribution, pipeline region-matched arm", spec["gain_built"]),
            ("rows retained, pipeline region-matched arm", spec["retained_built"]),
            ("share of the neg2 baseline excess over gc that is region mix",
             spec["region_share_of_excess"]),
            ("spearman(region-only AUROC, baseline rise, neg2 over gc)",
             spec["dose_baseline"]),
            ("spearman(region-only AUROC, contribution deficit, neg2 minus gc)",
             spec["dose_gain"])):
        v = must(label)
        if v is not None:
            near(label, v, key)

    # THE TWO GATES THAT CARRY THE DEFENCE. The composition-matched arms must be EXACTLY
    # uninformative on region -- 0.5 by construction, so anything else means a matcher
    # silently stopped matching region. And the ordering must survive removing region:
    # if it ever stops surviving, the bias-aware result is a region artefact and the
    # Discussion's mechanism paragraph is wrong rather than merely narrower.
    for arm in ("gc", "dn"):
        v = must(f"region-only AUROC, {arm} arm")
        if v is not None:
            record(abs(v - 0.5) <= spec["matched_arms_exact_tol"],
                   f"{arm} arm: region carries nothing, exactly, as its matcher requires",
                   f"{v:.6f}", "0.500000")
    # The rebuilt arm is held to the same standard: a stratified draw that lands anywhere but
    # 0.5 did not stratify.
    vb = must("region-only AUROC, pipeline region-matched arm")
    if vb is not None:
        record(abs(vb - 0.5) <= spec["matched_arms_exact_tol"],
               "the rebuilt arm matches region exactly, as a stratified draw must",
               f"{vb:.6f}", "0.500000")
    for label in ("region-matched neg2 still has the highest baseline",
                  "region-matched neg2 still has the lowest contribution"):
        v = must(label)
        if v is not None:
            record(int(v) == 1, label, int(v), 1)


def verify_peak_thresholds(T, g):
    """Were the positives filtered? ENCODE did it upstream; an earlier draft denied it."""
    print("\npeak thresholds  (ENCODE's released peaks are already filtered)")
    d = T.get("peak_thresholds.csv")
    if d is None:
        return record(False, "peak_thresholds.csv present", "MISSING",
                      "run scripts/peak_thresholds.py")
    spec = g["peak_thresholds"]
    q = d.set_index("check")

    def must(k):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, "value"])

    for label, key in (
            ("peak files audited", spec["n_files"]),
            ("peaks in the panel", spec["n_peaks"]),
            ("minimum log2 fold-enrichment over the panel", spec["min_log2fc"]),
            ("minimum -log10 p over the panel", spec["min_neglog10p"]),
            ("fraction of peaks carrying ENCODE's IDR label", spec["idr_fraction"])):
        v = must(label)
        if v is not None:
            near(label, v, key)

    # THE CLAIM THE LIMITATIONS NOW MAKES. Not "few peaks are weak" but "none is", so this
    # is a floor at the published threshold and a hard zero on violations.
    fc = must("minimum log2 fold-enrichment over the panel")
    if fc is not None:
        at_least("every peak clears ENCODE's fold-enrichment threshold", fc,
                 spec["fold_enrichment_floor"])
    pv = must("minimum -log10 p over the panel")
    if pv is not None:
        at_least("every peak clears ENCODE's significance threshold", pv,
                 spec["significance_floor"])
    bad = must("files below either threshold")
    if bad is not None:
        record(int(bad) == 0, "no peak file falls below either threshold", int(bad), 0)


def verify_match_quality(T, g):
    """What the matchers achieved. The paper argued from the specification for months."""
    print("\nmatch quality  (what the matchers DID, not what they were asked for)")
    d = T.get("match_quality.csv")
    if d is None:
        return record(False, "match_quality.csv present", "MISSING",
                      "run scripts/match_quality.py")
    spec = g["match_quality"]
    q = d.set_index("check")

    def must(k):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, "value"])

    n = q["n"].iloc[0] if "n" in q.columns else None
    if n is not None:
        record(int(n) == spec["n_datasets"], "datasets", int(n), spec["n_datasets"])
    for label, key in (
            ("gc_gap_median, gc arm", spec["gc_gap_median_gc"]),
            ("gc_gap_max, gc arm", spec["gc_gap_max_gc"]),
            ("gc_gap_median, dn arm", spec["gc_gap_median_dn"]),
            ("gc_gap_median, neg2 arm", spec["gc_gap_median_neg2"]),
            ("dinuc_l1_median, gc arm", spec["dinuc_l1_median_gc"]),
            ("dinuc_l1_median, dn arm", spec["dinuc_l1_median_dn"]),
            ("dinuc_l1_median, neg2 arm", spec["dinuc_l1_median_neg2"]),
            ("fraction of pairs within |dGC| 0.05, gc arm", spec["within_nominal_gc"]),
            ("dinucleotide L1 improvement factor, dn vs gc arm",
             spec["l1_improvement_factor"])):
        v = must(label)
        if v is not None:
            near(label, v, key)

    # THE FALLBACK CAP. The GC matcher accepts a best-seen candidate up to 3x its nominal
    # tolerance. If a pair ever exceeds that, the matcher is not the one described in Methods.
    mx = must("gc_gap_max, gc arm")
    if mx is not None:
        at_most("every GC pair is inside the matcher's real acceptance bound, which is the 3x "
                "fallback cap and NOT the nominal 0.05", mx, spec["gc_fallback_cap"])

    # THE RELAXATION LADDER MUST BE PRESENT, not merely correct. Reporting the nominal
    # tolerance alone overstates the match for 5.2% of pairs.
    if spec["ladder_must_be_reported"]:
        for rung in ("0.05", "0.10", "0.15"):
            record(f"fraction of pairs within |dGC| {rung}, gc arm" in q.index,
                   f"the GC relaxation ladder is reported at {rung}, so the nominal tolerance "
                   f"cannot be quoted alone", f"rung {rung}", "present")

    # AND THE DESIGN ARGUMENT MUST HOLD ON THE DATA. "1 of 15 vs 15 of 15 degrees of freedom"
    # is the reason the headline contrast's SIGN is design-implied; if the dinucleotide arm
    # were not materially better matched in fact, that argument would rest on nothing.
    imp = must("dinucleotide L1 improvement factor, dn vs gc arm")
    if imp is not None:
        record(imp >= spec["min_l1_improvement"],
               "the dinucleotide arm really is materially better composition-matched, so the "
               "1-of-15 vs 15-of-15 argument holds on the realisation and not only the design",
               f"{imp:.2f}x", f">= {spec['min_l1_improvement']}x")


def verify_score_scale(T, g):
    """R1t: the deep rungs entered on the probability scale, the k-mer on the logit scale."""
    print("\nR1t the score scale  (were the ladder's rungs measured the same way?)")
    d = T.get("score_scale_check.csv")
    if d is None:
        return record(False, "score_scale_check.csv present", "MISSING",
                      "run scripts/score_scale_check.py")
    spec = g["r1t_score_scale"]
    q = d.set_index("check")

    def must(k):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, "value"])

    n = q["n"].iloc[0] if "n" in q.columns else None
    if n is not None:
        record(int(n) == spec["n_datasets"], "datasets", int(n), spec["n_datasets"])
    rep = must("max |probability-scale gain - published gain|")
    if rep is not None:
        at_most("the probability-scale refit reproduces the PUBLISHED deep gain, so the "
                "deltas below are the scale and not a different pipeline",
                rep, spec["max_published_reproduction"])

    deltas = {}
    for label, key in (
            ("cnn logit minus probability scale, gc arm", spec["delta_cnn_gc"]),
            ("splicebert logit minus probability scale, gc arm", spec["delta_splicebert_gc"]),
            ("cnn logit minus probability scale, dn arm", spec["delta_cnn_dn"]),
            ("splicebert logit minus probability scale, dn arm", spec["delta_splicebert_dn"])):
        v = must(label)
        if v is not None:
            near(label, v, key)
            deltas[label] = v

    # THE BIAS IS REAL AND ONE-DIRECTIONAL, so the published deep gains are CONSERVATIVE.
    # Worth asserting rather than describing: if it ever changed sign, the published numbers
    # would be overstating the deep models and the section would have to be rewritten.
    if spec["scale_bias_must_be_one_directional"] and len(deltas) == 4:
        record(all(v > 0 for v in deltas.values()),
               "the probability scale UNDERSTATES every deep gain, so the published numbers "
               "are conservative for the deep models",
               f"min {min(deltas.values()):+.4f}", "all > 0")
        at_most("and the bias is small enough that the ladder's absolute values are not "
                "misstated", max(deltas.values()), spec["max_scale_bias"])

    # AND THE CLAIM IS INVARIANT. This is what makes it a disclosed sensitivity rather than a
    # correction: the contrast is the paper's quantity and it does not move.
    for label, key, pubkey in (
            ("cnn R1 contrast, logit scale", spec["contrast_cnn_logit"], "cnn"),
            ("splicebert R1 contrast, logit scale", spec["contrast_splicebert_logit"],
             "splicebert")):
        v = must(label)
        pv = must(f"{pubkey} R1 contrast, probability scale")
        if v is not None:
            near(label, v, key)
        if v is not None and pv is not None:
            at_most(f"the R1 contrast for the {pubkey} is INVARIANT to the score scale, so "
                    f"the paper's claim does not depend on the defect", abs(v - pv),
                    spec["max_contrast_shift"])


def verify_multiplier_variance(T, g):
    """R1s: the multiplier is mostly the protein's, but only against each factor's own null."""
    print("\nR1s whose property is the multiplier?  (and is 79 levels doing the work?)")
    d = T.get("multiplier_variance.csv")
    if d is None:
        return record(False, "multiplier_variance.csv present", "MISSING",
                      "run scripts/multiplier_variance.py")
    spec = g["r1s_multiplier_variance"]
    q = d.set_index("check")

    def must(k):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, "value"])

    n = q["n"].iloc[0] if "n" in q.columns else None
    if n is not None:
        record(int(n) == spec["n_cells"], "(dataset x model) cells, both arms positive",
               int(n), spec["n_cells"])
    for label, key in (
            ("share of log-multiplier variance, protein", spec["share_protein"]),
            ("share of log-multiplier variance, model", spec["share_model"]),
            ("share of log-multiplier variance, cell", spec["share_cell"]),
            ("permutation null share, protein", spec["null_protein"]),
            ("excess over the permutation null, protein", spec["excess_protein"]),
            ("panel multiplier, exp(mean log)", spec["panel_multiplier"])):
        v = must(label)
        if v is not None:
            near(label, v, key)

    # THE NULL MUST BE PRESENT, not merely correct. A decomposition reported without it is the
    # claim that was retracted here: 79 levels beat 3 levels on degrees of freedom alone.
    if spec["null_must_be_reported"]:
        for k in ("protein", "model", "cell"):
            record(f"permutation null share, {k}" in q.index,
                   f"the permutation null is reported for {k}, so its share cannot be read "
                   f"as a finding on its own", f"permutation null share, {k}", "present")

    # THE BLOCK-PRESERVING NULL AND THE DIRECT TEST. The wholesale null overstated the excess
    # five-fold; both corrections are gated so the generous number cannot be quoted alone.
    for label, key in (
            ("share of log-multiplier variance, dataset", spec["share_dataset"]),
            ("block-preserving null share, protein", spec["null_protein_block"]),
            ("excess over the block-preserving null, protein", spec["excess_protein_block"]),
            ("cross-cell-line correlation of the log multiplier", spec["cross_cell_r"]),
            ("cross-cell-line spearman of the log multiplier", spec["cross_cell_rho"])):
        v = must(label)
        if v is not None:
            near(label, v, key)
    k = "cross-cell-line correlation of the log multiplier"
    if k in q.index and "n" in d.columns:
        record(int(q.loc[k, "n"]) == spec["cross_cell_n"]["value"],
               "the cross-cell-line test uses every protein x model pair available from the "
               "proteins assayed in both lines", int(q.loc[k, "n"]),
               spec["cross_cell_n"]["value"])
    if spec["block_null_must_be_reported"]:
        record("excess over the block-preserving null, protein" in q.index,
               "the block-preserving excess is reported, so the five-fold-too-generous "
               "wholesale excess cannot be quoted on its own",
               "excess over the block-preserving null", "present")

    # B13. THE DIRECT TEST'S INTERVAL. The collapsed correlation is quoted as the paper's
    # cross-cell-line evidence, and at n = 15 a point estimate with a p and no interval invites
    # the reader to treat 0.598 as the finding. Both intervals are gated, and so is the fact
    # that the bootstrap one CONTAINS ZERO -- that is what licenses the text saying the result
    # is consistent with reproducibility rather than evidence for it.
    kc = "cross-cell-line correlation, collapsed over models"
    kb = "cross-cell-line correlation, collapsed, bootstrap CI"
    v = must(kc)
    if v is not None:
        near(kc, v, spec["cross_cell_collapsed_r"])
    if kc in q.index and "ci_low" in d.columns:
        near("cross-cell-line Fisher z lower bound",
             float(q.loc[kc, "ci_low"]), spec["cross_cell_fisher_low"])
        near("cross-cell-line Fisher z upper bound",
             float(q.loc[kc, "ci_high"]), spec["cross_cell_fisher_high"])
    if kb in q.index and "ci_low" in d.columns:
        lo, hi = float(q.loc[kb, "ci_low"]), float(q.loc[kb, "ci_high"])
        near("cross-cell-line bootstrap lower bound", lo, spec["cross_cell_boot_low"])
        near("cross-cell-line bootstrap upper bound", hi, spec["cross_cell_boot_high"])
        if spec["cross_cell_boot_must_span_zero"]:
            record(lo < 0 < hi,
                   "the cross-cell-line bootstrap interval CONTAINS ZERO at n=15, so the "
                   "direct test is consistent with reproducibility and is not evidence for it",
                   f"[{lo:+.3f}, {hi:+.3f}]", "spans zero")
        pe = must("excess over the permutation null, protein")
        pb = must("excess over the block-preserving null, protein")
        if pe is not None and pb is not None:
            record(pb < pe, "and it is the SMALLER of the two, which is why the wholesale "
                            "null was the wrong one to lead with",
                   f"block {pb:.3f} vs wholesale {pe:.3f}", "block < wholesale")

    if spec["protein_must_beat_its_own_null"]:
        s = must("share of log-multiplier variance, protein")
        nl = must("block-preserving null share, protein")
        if s is not None and nl is not None:
            record(s > nl, "protein's share beats the share it gets when protein labels are "
                           "permuted between datasets, so RBP identity carries information "
                           "beyond the individual experiment",
                   f"{s:.3f}", f"> {nl:.3f}")


def verify_horlacher(T, g):
    """R1p: the range replicates on an independent benchmark; R1n's strong form does not."""
    print("\nR1p Horlacher's benchmark  (does any of this travel off our own windows?)")
    d = T.get("horlacher_arm.csv")
    if d is None:
        return record(False, "horlacher_arm.csv present", "MISSING",
                      "run scripts/horlacher_arm.py")
    spec = g["r1p_horlacher"]
    q = d.set_index("check")

    def must(k):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, "value"])

    n = q["n"].iloc[0] if "n" in q.columns else None
    if n is not None:
        record(int(n) == spec["n_datasets"], "datasets from their release", int(n),
               spec["n_datasets"])
    for k, gk in (("composition alone, negative-1", "composition_n1"),
                  ("composition alone, negative-2", "composition_n2"),
                  ("nested contribution, negative-1", "gain_n1"),
                  ("nested contribution, negative-2", "gain_n2"),
                  ("CONTRAST, negative-2 minus negative-1", "contrast_n2_minus_n1"),
                  ("fold range across their two negative sets", "fold_range_their_data"),
                  ("within-dataset spearman(delta baseline, delta gain), their data",
                   "within_dataset_spearman")):
        v = must(k)
        if v is not None:
            near(k, v, spec[gk])

    fr = must("fold range across their two negative sets")
    if fr is not None:
        at_least("the RANGE replicates on an independent benchmark, with their positives, "
                 "their peaks and their folds", fr, spec["min_fold_range_their_data"])
    r = must("within-dataset spearman(delta baseline, delta gain), their data")
    if r is not None and spec["gradient_must_replicate_in_sign"]:
        record(r < 0, "the baseline gradient replicates in SIGN on their data",
               f"{r:+.3f}", "< 0")

    hi = must("n2 minus n1 gain, negative-2 baseline HIGHER")
    lo = must("n2 minus n1 gain, negative-2 baseline LOWER")
    if hi is not None:
        near("gain difference where negative-2 raises the baseline", hi,
             spec["gain_diff_n2_baseline_higher"])
    if lo is not None:
        near("gain difference where negative-2 lowers the baseline", lo,
             spec["gain_diff_n2_baseline_lower"])
    if hi is not None and lo is not None and spec["sign_must_not_reverse_on_their_data"]:
        record(hi < 0 and lo < 0,
               "THE LIMITATION: on their benchmark the sign does NOT reverse, so a "
               "protocol-specific residual remains and R1n's strong form is OURS only",
               f"{hi:+.4f}, {lo:+.4f}", "both negative")
        record(lo > hi, "the deficit still shrinks in the direction the baseline predicts",
               f"{lo:+.4f} > {hi:+.4f}", "gradient present")


def verify_recommendation(T, g):
    """R1q: the deliverable is tested. Direction is consistent; significance is one pair."""
    print("\nR1q recommendation  (does reporting the baseline actually help a reader?)")
    d = T.get("recommendation_works.csv")
    if d is None:
        return record(False, "recommendation_works.csv present", "MISSING",
                      "run scripts/recommendation_works.py")
    spec = g["r1q_recommendation"]
    q = d.set_index("check")

    def must(k, col="value"):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, col])

    for label, key in (
            ("rank agreement, raw, gc vs dn", spec["rank_raw_gc_dn"]),
            ("rank agreement, headroom, gc vs dn", spec["rank_headroom_gc_dn"]),
            ("rank agreement gain, gc vs dn", spec["rank_gain_gc_dn"]),
            ("scale-free disagreement, raw, gc vs dn", spec["disagreement_raw_gc_dn"]),
            ("scale-free disagreement, headroom, gc vs dn", spec["disagreement_headroom_gc_dn"]),
            ("scale-free disagreement, raw, dn vs neg2", spec["disagreement_raw_dn_neg2"]),
            ("scale-free disagreement, headroom, dn vs neg2",
             spec["disagreement_headroom_dn_neg2"])):
        v = must(label)
        if v is not None:
            near(label, v, key)

    npairs = must("number of protocol pairs")
    if npairs is not None:
        record(int(npairs) == spec["n_pairs"], "protocol pairs compared", int(npairs),
               spec["n_pairs"])
    ri = must("protocol pairs where rank agreement improves")
    di = must("protocol pairs where disagreement shrinks")
    if ri is not None:
        record(int(ri) == spec["pairs_rank_improves"],
               "pairs where rank agreement improves", int(ri), spec["pairs_rank_improves"])
    if di is not None:
        record(int(di) == spec["pairs_disagreement_shrinks"],
               "pairs where disagreement shrinks", int(di), spec["pairs_disagreement_shrinks"])
    if ri is not None and di is not None and spec["recommendation_must_help_on_every_pair"]:
        record(int(ri) == int(npairs or 3) and int(di) == int(npairs or 3),
               "THE DELIVERABLE IS TESTED: reporting the baseline helps on every protocol "
               "pair, on both measures", f"rank {int(ri)}/3, disagreement {int(di)}/3",
               "3/3 and 3/3")

    # THE LIMIT, asserted so the wording cannot drift stronger than the evidence.
    sig = 0
    for a, b in (("gc", "dn"), ("gc", "neg2"), ("dn", "neg2")):
        lo = must(f"rank agreement gain, {a} vs {b}", "ci_low")
        if lo is not None and lo > 0:
            sig += 1
    at_least("protocol pairs whose rank improvement is individually significant", sig,
             spec["min_significant_rank_pairs"])
    record(sig < 3, "NOT every pair's rank improvement is individually significant, so the "
                    "claim is direction-consistent rather than uniformly significant",
           f"{sig}/3", "< 3")


def verify_k_sweep(T, g):
    """R1 rebuilt from sequence, and shown not to depend on the k-mer size."""
    print("\nR1e  k sweep and rebuild-from-sequence")
    d = T.get("k_sweep.csv")
    if d is None:
        return record(False, "k_sweep.csv present", "MISSING", "run scripts/k_sweep.py")
    spec = g["r1_k_sweep"]
    q = d.set_index("check")

    def must(k):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, "value"])

    n = float(q.loc["contrast, k=4", "n"]) if "contrast, k=4" in q.index else None
    if n is not None:
        record(n == spec["n_datasets"]["value"], "datasets rebuilt from sequence", n,
               spec["n_datasets"]["value"])
    # Keys spelled out rather than built with an f-string. tests/unit/test_golden_keys_are_read
    # reads this file as TEXT, so `spec[f"contrast_k{k}"]` is invisible to it and four gated
    # keys registered as unread. A check the key-coverage test cannot see is a check that can
    # silently disappear later.
    for k, gk in ((3, "contrast_k3"), (4, "contrast_k4"),
                  (5, "contrast_k5"), (6, "contrast_k6")):
        v = must(f"contrast, k={k}")
        if v is not None:
            near(f"contrast, k={k}", v, spec[gk])

    diff = must("absolute difference")
    if diff is not None:
        at_most("REBUILT from sequence matches the committed contrast", diff,
                spec["max_rebuild_diff"])
    lo = must("smallest contrast across k=3..6")
    if lo is not None:
        at_least("contrast survives every k", lo, spec["min_contrast_any_k"])
    cnt = must("datasets positive at EVERY k")
    if cnt is not None:
        at_least("datasets positive at every k", cnt, spec["min_datasets_positive_all_k"])
    d54 = must("k=5 minus k=4")
    if d54 is not None:
        at_most("k=5 and k=4 agree (the '5-mer' error was inconsequential)", abs(d54),
                spec["max_k5_minus_k4"])


def verify_r1_robustness(T, g):
    """Replication across cell lines, and the efficiency the protocol buys.

    These are the paper's answer to its own biggest concession, that the sign of the contrast
    is design-implied. Neither is implied by the design.
    """
    print("\nR1d  replication and efficiency")
    d = T.get("r1_robustness.csv")
    if d is None:
        return record(False, "r1_robustness.csv present", "MISSING",
                      "run scripts/r1_robustness.py")
    spec = g["r1_robustness"]
    q = d.set_index("check")

    def must(k, col="value"):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, col])

    for k, gk in (("proteins assayed in both cell lines", "n_proteins_both_lines"),
                  ("REPLICATION of the contrast across cell lines", "replication_r"),
                  ("replication of the GC-arm gain alone", "replication_gc_alone"),
                  ("replication of the dinuc-arm gain alone", "replication_dn_alone"),
                  ("contrast r minus the better arm's r", "replication_ordering_gap"),
                  ("replication, PARTIALLING OUT total nested gain", "replication_partial"),
                  ("correlation of the contrast with total gain", "contrast_total_gain_corr"),
                  ("between-dataset sd of the contrast", "between_dataset_sd"),
                  ("p for the partialled replication", "partial_p"),
                  ("EFFICIENCY GAIN, z ratio", "efficiency_z_ratio"),
                  ("relative sample size, ratio of means", "relative_sample_size"),
                  ("relative sample size, median dataset", "relative_sample_median"),
                  ("relative sample size, mean over datasets", "relative_sample_mean"),
                  ("datasets needing MORE windows under dinuc matching",
                   "datasets_needing_more_windows"),
                  ("datasets where composition BEATS the model, GC arm",
                   "composition_beats_model_gc"),
                  ("datasets where composition BEATS the model, dinuc arm",
                   "composition_beats_model_dn")):
        v = must(k)
        if v is not None:
            near(k, v, spec[gk])

    lo = must("REPLICATION of the contrast across cell lines", "ci_low")
    if lo is not None:
        at_least("replication interval is well clear of zero", lo,
                 spec["min_replication_ci_low"])

    # The ordering must remain unassertable.
    olo = must("contrast r minus the better arm's r", "ci_low")
    ohi = must("contrast r minus the better arm's r", "ci_high")
    if olo is not None and ohi is not None and spec["ordering_ci_must_include_zero"]:
        record(olo <= 0 <= ohi, "replication ordering is NOT established (interval straddles 0)",
               f"[{olo:+.3f}, {ohi:+.3f}]", "includes 0")

    # THE LIMITATION, ASSERTED. The partial correlation must NOT be establishable, or R1d is
    # claiming more than the data support and the text has to change deliberately.
    plo = must("replication, PARTIALLING OUT total nested gain", "ci_low")
    phi = must("replication, PARTIALLING OUT total nested gain", "ci_high")
    if plo is not None and phi is not None and spec["partial_ci_must_include_zero"]:
        record(plo <= 0 <= phi,
               "replication net of total signal is NOT established (straddles 0)",
               f"[{plo:+.3f}, {phi:+.3f}]", "includes 0")

    z = must("EFFICIENCY GAIN, z ratio")
    if z is not None:
        at_least("proper matching is more statistically efficient", z,
                 spec["min_efficiency_ratio"])


def verify_strand_asymmetry(T, g):
    """Which arm carries the strand cue. Its direction bounds the artifact's effect."""
    print("\nR1  strand asymmetry  (which arm carries the cue?)")
    d = T.get("strand_asymmetry.csv")
    if d is None:
        return record(False, "strand_asymmetry.csv present", "MISSING",
                      "run scripts/strand_asymmetry.py")
    spec = g["r1_strand_asymmetry"]
    q = d.set_index("check")

    def must(k):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, "value"])

    for k, gk in (("frac_sense, GC arm", "frac_sense_gc"),
                  ("frac_sense, dinuc arm", "frac_sense_dn"),
                  ("asymmetry, dinuc minus GC", "asymmetry")):
        v = must(k)
        if v is not None:
            near(k, v, spec[gk])

    for k, gk in (("genes overlapping a sense-KEPT negative", "ngenes_kept"),
                  ("genes overlapping a DROPPED negative", "ngenes_dropped")):
        v = must(k)
        if v is not None:
            near(k, v, spec[gk])
    nk, nd = (must("genes overlapping a sense-KEPT negative"),
              must("genes overlapping a DROPPED negative"))
    if nk is not None and nd is not None and spec["kept_must_be_less_dense"]:
        record(nk < nd, "retention selects against multi-gene loci (justifies the stratum)",
               f"{nk:.3f} vs {nd:.3f}", "kept lower")

    asym = must("asymmetry, dinuc minus GC")
    if asym is not None and spec["asymmetry_must_be_positive"]:
        record(asym > 0, "cue is WEAKER in the GC arm, so the contrast is conservative",
               f"{asym:+.4f}", "> 0")
    cnt = must("datasets where dinuc is more sense")
    if cnt is not None:
        at_least("datasets where the dinuc arm is more sense", cnt,
                 spec["min_datasets_dinuc_more_sense"])


def verify_strand_placebo(T, g):
    """The pre-registered strand test: restriction against a REGION-MATCHED placebo.

    The criteria in golden.yaml were fixed before the experiment ran and must
    not be loosened afterwards. The one thing that DID change is a retraction: with an
    unstratified placebo the excess excluded zero and was reported as a real artifact; matched
    on region it does not, so only a bound is asserted now.
    """
    print("\nR1  strand placebo  (the pre-registered test)")
    d = T.get("strand_placebo.csv")
    if d is None:
        return record(False, "strand_placebo.csv present", "MISSING",
                      "run scripts/strand_placebo.py")
    spec = g["r1_strand_placebo"]
    q = d.set_index("check")

    def must(k, col="value"):
        if k not in q.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(q.loc[k, col])

    n = must("contrast, full data", "n")
    if n is not None:
        record(n == spec["n_datasets"]["value"], "datasets in the test", n,
               spec["n_datasets"]["value"])

    for k, gk in (("contrast, full data", "contrast_full"),
                  ("contrast, sense-only pairs", "contrast_sense_only"),
                  ("contrast, PLACEBO (same n, random)", "contrast_placebo"),
                  ("contrast, PLACEBO stratified on region x GC", "contrast_placebo_strat"),
                  ("change from restriction", "change_from_restriction"),
                  ("change from placebo", "change_from_placebo"),
                  ("strand excess, UNSTRATIFIED placebo", "strand_excess_unstratified"),
                  ("locus-mix component", "locus_mix"),
                  ("STRAND-SPECIFIC EXCESS (stratified)", "strand_excess"),
                  ("strand-CORRECTED contrast", "corrected_contrast")):
        v = must(k)
        if v is not None:
            near(k, v, spec[gk])

    ex = must("STRAND-SPECIFIC EXCESS (stratified)")
    if ex is not None:
        at_most("strand artifact is bounded and small", abs(ex),
                spec["max_abs_strand_excess"])

    # The locus mix must be SMALL. At 5 seeds it looked real (-0.0024, interval clear of zero)
    # and this gate demanded that it be. At 20 seeds it is +0.0004 and indistinguishable from
    # zero, so the earlier gap was largely seed noise and the assertion is now a bound.
    lm = must("locus-mix component")
    if lm is not None:
        at_most("locus mix is small once the placebo has enough seeds", abs(lm),
                spec["max_abs_locus_mix"])

    # THE ARTIFACT IS REAL. Withdrawn once on a 5-seed interval that touched zero; restored at
    # 20 seeds. Asserted so it cannot be withdrawn again without the numbers moving.
    elo = must("STRAND-SPECIFIC EXCESS (stratified)", "ci_low")
    ehi = must("STRAND-SPECIFIC EXCESS (stratified)", "ci_high")
    if elo is not None and ehi is not None and spec["excess_ci_must_exclude_zero"]:
        record(not (elo <= 0 <= ehi), "strand artifact is REAL (interval excludes zero)",
               f"[{elo:+.4f}, {ehi:+.4f}]", "excludes 0")

    # THE PRE-REGISTERED CRITERIA.
    clo = must("strand-CORRECTED contrast", "ci_low")
    cv = must("strand-CORRECTED contrast")
    if cv is not None and clo is not None and spec["corrected_ci_must_exclude_zero"]:
        record(cv > 0 and clo > 0, "PRE-REGISTERED: sign kept and CI excludes zero",
               f"{cv:+.4f} [{clo:+.4f}, ...]", "> 0")
    fr = must("fraction of the contrast surviving")
    if fr is not None:
        at_least("PRE-REGISTERED: fraction of the contrast surviving", fr,
                 spec["min_fraction_surviving"])


def verify_unconditional_refit(T, g):
    """The corrected attenuation analysis, and its calibrated null.

    Gated separately from r4_incremental_value because it CONTRADICTS it. The retracted
    version read near-zero attenuation as independence; with a proper null that same number
    means substantial sharing. Both sets of keys exist so the contradiction is visible in the
    config rather than resolved silently.

    Missing rows are failures, per verify_scale_check.
    """
    print("\nunconditional refit  (near-zero attenuation is NOT independence)")
    d = T.get("unconditional_refit.csv")
    if d is None:
        return record(False, "unconditional_refit.csv present", "MISSING",
                      "run scripts/unconditional_refit.py")
    spec = g["r4_unconditional_refit"]
    w = d[d.standardisation == "within_dataset"].set_index("check")
    allrows = d.set_index(["standardisation", "check"])

    def must(k):
        if k not in w.index:
            record(False, f"row present: {k}", "MISSING", "the row")
            return None
        return float(w.loc[k, "value"])

    for k, gk in (("coef, TRULY unconditional", "coef_unconditional_within"),
                  ("coef, conditional on phyloP", "coef_conditional_within"),
                  ("NULL attenuation at rho=0 (simulated)", "null_attenuation_simulated"),
                  ("NULL attenuation at rho=0 (analytic reference only)",
                   "null_attenuation_analytic"),
                  ("correlation implied by the observed attenuation", "implied_rho"),
                  ("spearman(|delta|, phyloP), MEASURED", "measured_rho")):
        val = must(k)
        if val is not None:
            near(k, val, spec[gk])

    # The raw attenuation is small and sign-unstable. Asserted as small so that a run where
    # it is LARGE is flagged: that would mean the estimator changed, not that the claim got
    # stronger.
    raw = must("attenuation fraction (identical rows)")
    if raw is not None:
        at_most("raw attenuation is near zero (and not the claim)", abs(raw),
                spec["max_abs_raw_attenuation"])

    # NO AGREEMENT ASSERTION between the two routes, deliberately. The closed form is
    # anti-conservative at this covariate strength, so gating agreement would gate a
    # known-invalid approximation against the correct simulation. Both are pinned separately.

    # THE CLAIM, required under BOTH standardisations because the raw attenuation flips sign
    # between them and the conclusion must not.
    for tag in ("pooled", "within_dataset"):
        k = (tag, "excess attenuation over the null")
        if k not in allrows.index:
            record(False, f"row present: {tag} excess attenuation", "MISSING", "the row")
            continue
        at_least(f"excess attenuation over the null, {tag}",
                 float(allrows.loc[k, "value"]), spec["min_excess_attenuation"])

    # ...and the sharing the attenuation implies must match the sharing measured directly.
    imp = must("correlation implied by the observed attenuation")
    mea = must("spearman(|delta|, phyloP), MEASURED")
    if imp is not None and mea is not None:
        at_most("implied and measured correlation agree", abs(imp - mea),
                spec["max_implied_measured_gap"])


def verify_strand_audit(T, g):
    """The negatives are ~45% antisense, and that must not be what the contrast measures.

    THIS FUNCTION EXISTS BECAUSE ITS ABSENCE WAS A FALSE STATEMENT. `golden.yaml` grew a
    `strand_audit` block with 9 keys that nothing read, while the write-up claimed it was gated. It
    was not: `grep strand scripts/verify.py` returned nothing. That is the 27th unread golden
    key in this project, created by the commit that wired up the first 26 -- the same bug
    class, one layer up, for the third time. The lesson that finally stuck is mechanical, not
    intellectual: `tests/unit/test_golden_keys_are_read.py` now fails the build if any key in
    golden.yaml is unreferenced, so this cannot recur by inspection failure again.
    """
    print("\nstrand audit  (the negatives are ~45% antisense; does it explain the contrast?)")
    r = T.get("strand_audit.csv")
    s = T.get("strand_audit_summary.csv")
    if r is None or s is None:
        return record(False, "strand audit tables present", "MISSING", "the tables")
    spec = g["strand_audit"]

    near("fraction of negatives on their own gene's strand", float(r.frac_sense.mean()),
         spec["frac_sense"])
    at_most("negatives outside any annotated gene", float(r.frac_no_gene.mean()),
            spec["frac_no_gene_max"])

    # SIX CHECKS WERE DELETED HERE, DELIBERATELY. They asserted that the strand artifact does
    # not predict the composition-SHARE contrast across antisense-rich and antisense-poor
    # halves. That contrast is retracted -- it is an algebraic identity -- so those six checks
    # were defending a claim the paper no longer makes while still counting toward the passing
    # total. `s` (strand_audit_summary.csv) is still required to exist above, because its
    # absence would mean the audit did not run, but nothing in it is asserted any more.
    del s


def verify_recompute(T, g):
    """The only check here that proves a number rather than reproducing it.

    Every other assertion in this file reads a table the analysis pass wrote, so it detects
    drift and not error: corrupting an input and re-running verify.py produces zero failures.
    This one reads the per-example model scores committed under data/evidence/ and recomputes
    the published AUROC from them, for both deep architectures. Zeroing the scores makes it
    report max|diff| ~0.45 and fail, so it has power.

    Kept deliberately duplicative of scripts/recompute.py: that script is what a cloner runs
    in three seconds with no credentials, and this is the same arithmetic inside the gate.
    """
    print("\nrecompute  (published AUROCs, rebuilt from committed per-example scores)")
    spec = g.get("recompute")
    if spec is None:
        return
    four = T.get("matched_four_models.csv")
    ev = ROOT / "data" / "evidence" / "scores"
    if four is None:
        return record(False, "matched_four_models.csv present", "MISSING", "the table")
    if not ev.exists():
        return record(False, "committed per-example scores present", "MISSING", str(ev))

    import glob as _glob

    from sklearn.metrics import roc_auc_score
    pub = four.set_index("dataset")
    for model, col in (("splicebert", "splicebert"), ("cnn", "cnn")):
        if col not in pub.columns:
            continue
        diffs = []
        for ds, want in pub[col].items():
            protein, cell = ds.split(":")
            fs = sorted(_glob.glob(str(ev / cell / protein / model / "fold*" / "scores.tsv.gz")))
            if len(fs) != 5:
                continue
            d = pd.concat([pd.read_csv(f, sep="\t") for f in fs], ignore_index=True)
            if d.label.nunique() < 2:
                continue
            diffs.append(abs(roc_auc_score(d.label, d.score) - float(want)))
        if not diffs:
            record(False, f"{model} recomputed from evidence", "NO DATASETS", ">= 1")
            continue
        at_least(f"{model}: datasets recomputed from evidence", len(diffs),
                 spec["min_datasets"])
        at_most(f"{model}: max |recomputed - published|", max(diffs), spec["max_abs_diff"])


def verify_cache_evidence(T, g):
    """The per-dataset tables the --from-cache paths read are EVIDENCE, and were ungated.

    THIS EXISTS BECAUSE AN ATTACK GOT THROUGH, AGAIN. run.sh regenerates five summaries with
    --from-cache, which reads a committed *_per_dataset.csv rather than redoing the refits.
    `grep per_dataset scripts/verify.py` returned nothing, so those tables were load-bearing and
    unasserted. Zeroing every per-arm gain column in k_sweep_per_dataset.csv and every AUROC
    column in strand_placebo_per_dataset.csv, then rebuilding, reproduced both summaries
    BIT-FOR-BIT and passed. R1e's claim to have rebuilt the headline from raw sequence therefore
    reproduced from a table in which all 188 per-arm gains were zero.

    The fix is to assert the arithmetic that links the evidence columns to the summary columns.
    A derived column that is never checked against what derives it is not evidence.
    """
    print("\ncache evidence  (the per-dataset tables --from-cache reads)")
    tol = g["integrity"]["max_cache_arithmetic_diff"]

    k = T.get("k_sweep_per_dataset.csv")
    if k is None:
        record(False, "k_sweep_per_dataset.csv present", "MISSING", "the table")
    else:
        for kk in (3, 4, 5, 6):
            c, gc, dn = f"contrast_k{kk}", f"gain_gc_k{kk}", f"gain_dn_k{kk}"
            if not {c, gc, dn} <= set(k.columns):
                record(False, f"k_sweep columns for k={kk}", "MISSING", "present")
                continue
            at_most(f"k={kk}: contrast equals dinuc gain minus GC gain",
                    float((k[c] - (k[dn] - k[gc])).abs().max()), tol)
        # ...and the k=4 gains must be the published ones, or the "rebuilt from sequence"
        # claim is a column read.
        gcp, dnp = T.get("rehearsal_binding_gc.csv"), T.get("rehearsal_binding_dinuc.csv")
        if gcp is not None and dnp is not None and "gain_gc_k4" in k.columns:
            for arm, pub in (("gc", gcp), ("dn", dnp)):
                j = k.merge(pub, on="dataset", how="inner")
                at_most(f"k=4 {arm} gain matches the published delta_auroc",
                        float((j[f"gain_{arm}_k4"] - j["delta_auroc"]).abs().max()),
                        g["integrity"]["max_rebuild_vs_published_diff"])

    # STRAND PLACEBO: derive the summary FROM the evidence columns and check it against what
    # the summary table reports. Checking a derived column against its own components would
    # pass on a table where both are zero; recomputing the reported means from the per-arm
    # AUROCs is what makes zeroing the evidence fail.
    sp = T.get("strand_placebo_per_dataset.csv")
    summ = T.get("strand_placebo.csv")
    if sp is None or summ is None:
        record(False, "strand_placebo per-dataset and summary present", "MISSING", "both")
    else:
        q = summ.set_index("check")
        pairs = (("contrast, full data", "full_dn", "full_gc"),
                 ("contrast, sense-only pairs", "sense_dn", "sense_gc"),
                 ("contrast, PLACEBO (same n, random)", "placebo_dn", "placebo_gc"),
                 ("contrast, PLACEBO stratified on region x GC",
                  "placebo_strat_dn", "placebo_strat_gc"))
        for row, dn, gc in pairs:
            if not {dn, gc} <= set(sp.columns) or row not in q.index:
                record(False, f"evidence for '{row}'", "MISSING", f"{dn}, {gc}")
                continue
            derived = float((sp[dn] - sp[gc]).mean())
            at_most(f"summary '{row[:34]}' derives from its per-arm evidence",
                    abs(derived - float(q.loc[row, "value"])), tol)
        # ...and that evidence must itself be the published numbers, not zeros.
        if "full_gc" in sp.columns and gcp is not None:
            j = sp.merge(gcp, on="dataset", how="inner")
            at_most("strand placebo's full GC gain matches the published delta_auroc",
                    float((j["full_gc"] - j["delta_auroc"]).abs().max()),
                    g["integrity"]["max_rebuild_vs_published_diff"])


def verify_cross_tables(T, g):
    """The two tables holding the SAME numbers must agree.

    THIS EXISTS BECAUSE AN ATTACK GOT THROUGH. verify.py gates R1 on cost_of_matching.csv and
    never opened rehearsal_binding_gc.csv or rehearsal_binding_dinuc.csv, which are where those
    numbers come from and what the manuscript's component means and scale_check.py both read.
    Nothing asserted the two agreed. Permuting rehearsal_binding_dinuc.csv against its dataset
    labels therefore passed 166/166 while turning "larger in 88/94" into 67/94, "94/94 fall"
    into 80/94, and the Wilcoxon p from 3.8e-17 into 1.5e-12 -- the last of which violates this
    file's own wilcoxon_p_max of 1e-15, in a check that never ran on the corrupted table.

    A duplicated number is only as trustworthy as the assertion that the copies match.
    """
    print("\ncross-table consistency")
    cm = T.get("cost_of_matching.csv")
    gc = T.get("rehearsal_binding_gc.csv")
    dn = T.get("rehearsal_binding_dinuc.csv")
    if cm is None or gc is None or dn is None:
        return record(False, "both rehearsal arms and cost_of_matching present", "MISSING",
                      "all three tables")
    tol = g["integrity"]["max_cross_table_diff"]
    for arm, src in (("gc", gc), ("dn", dn)):
        j = cm.merge(src, on="dataset", how="inner", suffixes=("", "_src"))
        record(len(j) == len(cm), f"{arm}: every cost_of_matching row present in the arm",
               len(j), len(cm))
        for cm_col, src_col in ((f"auroc_{arm}", "auroc"),
                                (f"composition_auroc_{arm}", "composition_auroc"),
                                (f"delta_auroc_{arm}", "delta_auroc")):
            if cm_col not in j or src_col not in j:
                record(False, f"{arm}: column {cm_col}", "MISSING", "present")
                continue
            at_most(f"{arm}: {cm_col} matches the rehearsal table",
                    float((j[cm_col] - j[src_col]).abs().max()), tol)


def verify_integrity(T, g):
    print("\nintegrity")
    spec = g["integrity"]

    # EVERY NUMBER IN THE MANUSCRIPT MUST HAVE A SOURCE. scripts/audit_manuscript.py lists the
    # ones that appear in no committed table and in no golden key. The paper's primary contrast
    # was in exactly that state through six rounds of adversarial review, because a reviewer
    # reads a number rather than goes looking for it. Ratcheted to zero, over both decimal
    # values and bare integers, so any unsourced number at all fails the build.
    orph = T.get("manuscript_orphans.csv")
    if orph is None:
        record(False, "manuscript_orphans.csv present", "MISSING",
               "run scripts/audit_manuscript.py")
    else:
        at_most("manuscript numbers with no source in any table", len(orph),
                spec["max_manuscript_orphans"])

    # NaN fraction across every published table. This key existed and was never read.
    #
    # COLUMN-AWARE, and the distinction is not a loophole. One column legitimately carries
    # NaN as a value: auroc_block_prevalence is undefined for a dataset whose 1-Mb blocks
    # never contain both a pathogenic and a benign variant, which happens for the 8 datasets
    # holding between 1 and 6 pathogenic variants. NaN there means "not computable", and
    # writing 0.5 or dropping the row would both be worse -- one invents a result, the other
    # hides that the baseline could not be run.
    #
    # None of those 8 are in the powered stratum, so no reported number depends on them. The
    # exemption is listed per column rather than per table so that a NaN appearing anywhere
    # else still fails.
    # The three *_common columns are NaN on EXACTLY the same rows and for exactly the same
    # reason: they exist to score the model on the variants the block-prevalence baseline can
    # score, so where the baseline is undefined they are undefined too. Added when
    # the paired comparison was moved onto a common mask; the gate caught them immediately,
    # the coverage columns were first written.
    NAN_OK = {("variant_specificity.csv", "auroc_block_prevalence"): "undefined below ~10 "
              "pathogenic variants; excluded from every reported stratum",
              ("variant_specificity.csv", "n_common"): "defined only where the baseline is",
              ("variant_specificity.csv", "auroc_matched_common"): "defined only where the "
              "baseline is; this is the point of the column",
              ("variant_specificity.csv", "auroc_conservation_common"): "defined only where "
              "the baseline is"}
    worst, where = 0.0, ""
    for name in ("cost_of_matching.csv", "matched_four_models.csv", "locality_ism.csv",
                 "variant_specificity.csv", "variant_ladder_paired.csv"):
        d = T.get(name)
        if d is None:
            continue
        num = d.select_dtypes(include="number")
        cols = [c for c in num.columns if (name, c) not in NAN_OK]
        if not cols:
            continue
        frac = float(num[cols].isna().to_numpy().mean())
        if frac > worst:
            worst, where = frac, name
    at_most(f"NaN fraction, claim-bearing columns ({where or 'all clean'})",
            worst, spec["max_nan_fraction"])
    # And the exempted column is reported rather than ignored, so it cannot quietly grow.
    d = T.get("variant_specificity.csv")
    if d is not None and "auroc_block_prevalence" in d:
        at_most("uncomputable block-prevalence baselines, fraction",
                float(d.auroc_block_prevalence.isna().mean()), 0.15)
    for name in ("cost_of_matching.csv", "locality_ism.csv"):
        d = T.get(name)
        if d is None:
            continue
        if "dataset" in d:
            record(int(d.dataset.duplicated().sum()) == spec["duplicate_datasets_allowed"],
                   f"no duplicate datasets in {name}",
                   int(d.dataset.duplicated().sum()), spec["duplicate_datasets_allowed"])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--local", default=None, help="read tables from a local directory")
    p.add_argument("--golden", default=str(ROOT / "config" / "golden.yaml"))
    p.add_argument("--strict", action="store_true")
    a = p.parse_args()

    g = yaml.safe_load(Path(a.golden).read_text())
    T = Tables(a.local)

    print("=" * 78)
    print("VERIFY -- does the reproduction match the original science?")
    # golden.yaml's meta keys are reference_run/established. They were source_run/measured_on
    # when the file described the earlier study, and renaming them there left this line
    # reading keys that no longer exist -- a KeyError before a single check ran.
    m = g["meta"]
    print(f"golden: {m['reference_run']} established {m['established']}")
    print("=" * 78)

    for fn in (verify_r1, verify_scale_check, verify_r2, verify_r3, verify_r4_paired, verify_r4,
               verify_multidonor, verify_incremental_value, verify_unconditional_refit,
               verify_strand_contrast, verify_region, verify_deep_contrast, verify_protocol_identification, verify_expression_control, verify_cluster_intervals, verify_three_arm, verify_baseline_confounding, verify_scale_sweep, verify_protocol_or_baseline, verify_baseline_order, verify_baseline_order_models, verify_models_by_protocol, verify_three_arm_models, verify_transport, verify_fold_integrity, verify_region_asymmetry, verify_peak_thresholds, verify_design_effect, verify_standalone_auroc, verify_match_quality, verify_score_scale, verify_multiplier_variance, verify_horlacher, verify_recommendation, verify_k_sweep, verify_r1_robustness, verify_strand_asymmetry, verify_strand_placebo,
               verify_strand_audit, verify_recompute,
               verify_estimands, verify_nested_scale,
               verify_cache_evidence, verify_cross_tables, verify_integrity):
        try:
            fn(T, g)
        except Exception as e:                  # a broken check is a failure, not a crash
            record(False, f"{fn.__name__} raised", type(e).__name__, "no exception", str(e)[:80])

    # HOW MANY CHECKS RAN IS ITSELF A CHECK, and it is the only one that closes the whole
    # class of silent skips at once. Most gates in this file are still written
    # `if value is not None:`, so a table that exists but has had rows removed makes them
    # vanish rather than fail -- which is exactly how deleting five rows once produced
    # "106/106 passed". Rewriting 38 call sites would fix the instances; asserting the count
    # fixes the category, including gates nobody has thought to attack yet.
    #
    # This floor MUST be raised whenever checks are legitimately added, the same discipline
    # test_golden_keys_are_read.py enforces for golden keys. That is the intended cost.
    floor = g["integrity"]["min_domain_checks"]
    n_ran = len(checks)                          # domain checks only; this one is not counted
    record(n_ran >= floor, "number of domain checks that ran", n_ran, f">= {floor}",
           "" if n_ran >= floor else "gates were SKIPPED, not passed -- look for missing rows")

    bad = [c for c in checks if not c[0]]

    # RECORD WHAT RAN, so the manuscript's claim about coverage is auditable. The paper says
    # "696 numeric assertions"; that integer lived in prose and in no table, so
    # audit_manuscript.py could not source it and it was checked by hand. Written on every
    # run, pass or fail, because a claim of 696 against a run of 690 is the interesting case.
    try:
        pd.DataFrame([{"name": "assertions", "value": len(checks)},
                      {"name": "assertions_passed", "value": len(checks) - len(bad)},
                      {"name": "domain_checks", "value": n_ran}]).to_csv(
            ROOT / "results" / "tables" / "verify_summary.csv", index=False)
    except OSError:
        pass                                     # a read-only checkout still verifies

    print("\n" + "=" * 78)
    print(f"{len(checks) - len(bad)}/{len(checks)} checks passed")
    if bad:
        print("\nFAILED CLAIMS:")
        for _, claim, got, want, note in bad:
            print(f"  {claim}: got {got}, want {want}" + (f" ({note})" if note else ""))
        print("\nThe pipeline ran but the science did not reproduce. Do not write it up.")
        return 1
    print("Every golden number reproduced within tolerance.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
