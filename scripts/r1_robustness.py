"""R1d: does the contrast replicate, and does it buy anything? Two questions, committed tables.

WHY THESE TWO. The manuscript's standing concession is that the SIGN of the contrast is
implied by the design -- the composition baseline's 15 degrees of freedom are exactly what the
dinucleotide matcher controls -- so only the MAGNITUDE is informative. That concession is
correct and it is also the paper's biggest weakness, because a referee reads "only the
magnitude is informative" and asks what the magnitude is worth. These are the two answers that
do not require any new data.

  REPLICATION   Fifteen proteins were assayed in BOTH HepG2 and K562. Those are separate eCLIP
                experiments with separately drawn negatives, so a per-protein contrast measured
                in one line is an out-of-sample prediction of the other. If the magnitude is a
                reproducible property of the protein rather than noise around a design-implied
                sign, it must correlate across lines. Nothing about the design forces that: the
                sign is guaranteed in both lines independently, so a high correlation is
                information about magnitude and not about direction.

  EFFICIENCY    A benchmark builder does not act on an effect size, they act on whether the
                protocol lets them detect something. The nested gain rises 2.5x under
                dinucleotide matching, but so does its standard error, and only the ratio
                matters. z = gain / SE is that ratio, and because z grows as sqrt(n), a gain in
                z converts directly into how many labelled windows the same conclusion needs.
                This is the one number in the paper a practitioner can act on, and it is
                immune to the design-implied-sign objection entirely: the sign being
                predictable says nothing about how efficiently the effect is measured.

SE is recovered as (hi - lo) / (2 * 1.96). THE PROVENANCE MATTERS AND AN EARLIER VERSION OF
THIS DOCSTRING GOT IT WRONG: it called the source a bootstrap percentile interval, for which
that recovery would only be approximate. It is not one. `nested.gain_over_composition` writes
DeLong's diff +/- z*se (see src/rbp/eval/nested.py), which is symmetric by construction --
verified at max asymmetry 2.0e-16 across all 188 rows -- so the recovery is exact. In a paper
whose header repudiates one false provenance claim, another sitting under the efficiency
headline is exactly what a referee finds.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, wilcoxon

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
TABLES = ROOT / "results" / "tables"
N_BOOT = 2000
SEED = 0
Z95 = 1.959963985


def log(m):
    print(m, flush=True)


def main():
    cm = pd.read_csv(TABLES / "cost_of_matching.csv")
    cm["contrast"] = cm.delta_auroc_dn - cm.delta_auroc_gc
    rows = []

    def add(check, value, lo=np.nan, hi=np.nan, n=len(cm), note=""):
        rows.append({"check": check, "value": float(value), "ci_low": lo, "ci_high": hi,
                     "n": n, "note": note})

    # --- replication across cell lines ---------------------------------------------------
    w = cm.pivot_table(index="protein", columns="cell", values="contrast")
    w = w.dropna()
    n_rep = len(w)
    a, b = w.iloc[:, 0].to_numpy(), w.iloc[:, 1].to_numpy()
    r = pearsonr(a, b)
    rs = spearmanr(a, b)
    rng = np.random.default_rng(SEED)
    boots = [pearsonr(a[i], b[i]).statistic
             for i in (rng.integers(0, n_rep, n_rep) for _ in range(N_BOOT))
             if len(np.unique(a[i])) > 2]
    lo, hi = np.percentile(boots, [2.5, 97.5])
    add(f"proteins assayed in both cell lines", float(n_rep), n=n_rep,
        note=f"{' vs '.join(w.columns)}; separate experiments, separately drawn negatives")
    add("REPLICATION of the contrast across cell lines", r.statistic, lo, hi, n_rep,
        note=f"pearson, p={r.pvalue:.3g}; spearman {rs.statistic:+.3f}")
    add("mean |difference| between the two lines", float(np.abs(a - b).mean()), n=n_rep,
        note=f"against a between-dataset sd of {cm.contrast.std():.4f}")

    # The components, for comparison. The contrast replicates AT LEAST AS WELL as either, and
    # the paper says only that: a paired protein bootstrap gives r_contrast - max(r_arm) =
    # +0.113 [-0.082, +0.465], P(<=0) = 0.21 on n = 15, so the ordering is a point estimate and
    # not a finding. An earlier version asserted the ordering and gated it.
    for col, lab in (("delta_auroc_gc", "GC-arm gain"), ("delta_auroc_dn", "dinuc-arm gain")):
        p = cm.pivot_table(index="protein", columns="cell", values=col).dropna()
        rr = pearsonr(p.iloc[:, 0], p.iloc[:, 1])
        add(f"replication of the {lab} alone", rr.statistic, n=len(p),
            note=f"p={rr.pvalue:.3g}; compare with the contrast above")

    # THE ORDERING, TESTED RATHER THAN ASSERTED. An earlier version claimed the contrast
    # replicates BETTER than either arm and gated it. A paired protein bootstrap says the
    # ordering is a point estimate, not a finding, so the text now says "at least as well as".
    pg = cm.pivot_table(index="protein", columns="cell", values="delta_auroc_gc").dropna()
    pd_ = cm.pivot_table(index="protein", columns="cell", values="delta_auroc_dn").dropna()
    common = w.index.intersection(pg.index).intersection(pd_.index)
    rng2 = np.random.default_rng(SEED)
    dif = []
    for _ in range(N_BOOT):
        idx = common[rng2.integers(0, len(common), len(common))]
        try:
            f = lambda P: pearsonr(P.loc[idx].iloc[:, 0], P.loc[idx].iloc[:, 1]).statistic
            dif.append(f(w) - max(f(pg), f(pd_)))
        except Exception:
            continue
    dif = np.asarray(dif)
    dlo, dhi = np.percentile(dif, [2.5, 97.5])
    add("contrast r minus the better arm's r", float(dif.mean()), dlo, dhi, len(common),
        note=f"P(<=0) = {float((dif <= 0).mean()):.2f}; the ordering is NOT established")

    # --- statistical efficiency ----------------------------------------------------------
    gc = pd.read_csv(TABLES / "rehearsal_binding_gc.csv")
    dn = pd.read_csv(TABLES / "rehearsal_binding_dinuc.csv")
    m = gc.merge(dn, on="dataset", suffixes=("_gc", "_dn"))
    for arm in ("gc", "dn"):
        m[f"se_{arm}"] = (m[f"delta_ci_high_{arm}"] - m[f"delta_ci_low_{arm}"]) / (2 * Z95)
        m[f"z_{arm}"] = m[f"delta_auroc_{arm}"] / m[f"se_{arm}"].replace(0, np.nan)
    m = m.dropna(subset=["z_gc", "z_dn"])
    wz = wilcoxon(m.z_dn, m.z_gc)
    ratio = m.z_dn.mean() / m.z_gc.mean()
    add("mean z = gain / SE, GC arm", float(m.z_gc.mean()), n=len(m),
        note=f"median {m.z_gc.median():.2f}")
    add("mean z = gain / SE, dinuc arm", float(m.z_dn.mean()), n=len(m),
        note=f"median {m.z_dn.median():.2f}")
    add("EFFICIENCY GAIN, z ratio", float(ratio), n=len(m),
        note=f"higher in {int((m.z_dn > m.z_gc).sum())}/{len(m)}, paired Wilcoxon "
             f"p={wz.pvalue:.3g}")
    # z grows as sqrt(n), so the sample needed for a fixed z scales as 1/ratio^2. REPORTED
    # THREE WAYS, because the ratio of means flatters it: the mean of the per-dataset ratios
    # is close to 1, and some datasets are worse off. Quoting only the first would be a
    # ratio-of-means artifact.
    per = (m.z_gc / m.z_dn) ** 2
    add("relative sample size, ratio of means", float(1.0 / ratio ** 2), n=len(m),
        note="z ~ sqrt(n); the headline form, and the most favourable")
    add("relative sample size, median dataset", float(per.median()), n=len(m), note="")
    add("relative sample size, mean over datasets", float(per.mean()), n=len(m),
        note="close to 1: the advantage is concentrated, not universal")
    add("datasets needing MORE windows under dinuc matching", float((per > 1).sum()),
        n=len(m), note="the protocol is not better everywhere, and the paper says so")
    # And the blunt version a reader will remember.
    add("datasets where composition BEATS the model, GC arm",
        float((m.composition_auroc_gc >= m.auroc_gc).sum()), n=len(m), note="")
    add("datasets where composition BEATS the model, dinuc arm",
        float((m.composition_auroc_dn >= m.auroc_dn).sum()), n=len(m),
        note="the harder protocol makes the model look better on this measure, not worse")

    out = pd.DataFrame(rows)
    out.to_csv(TABLES / "r1_robustness.csv", index=False)
    for _, x in out.iterrows():
        ci = f" [{x.ci_low:+.3f}, {x.ci_high:+.3f}]" if pd.notna(x.ci_low) else ""
        log(f"  {x.check:52} {x.value:+.4f}{ci}   {x.note}")
    log(f"\n  wrote {TABLES / 'r1_robustness.csv'}")


if __name__ == "__main__":
    main()
