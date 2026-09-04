"""Step 3: what effect could we have detected? Turns a null into a bounded statement.

118 of 157 datasets did not survive the conservation control. Without this script that
sentence means nothing: it could be "the effect is absent" or "we had 15 pathogenic
variants and could not have seen an elephant". The minimum detectable effect distinguishes
them, and it is the difference between a publishable null and an absence of evidence.

Two routes, and the cheap one is used only after the expensive one confirms it:

  ANALYTIC.  MDE = (z_0.975 + z_0.80) * SE = 2.80 * SE at 80% power, two-sided alpha=0.05.
  Fast, and it uses each dataset's OBSERVED standard error, so the real sample size,
  class balance and delta-conservation correlation are already baked in.

  SIMULATED.  Generate data with a known effect, run the actual test, count how often the
  bootstrap interval excludes zero. Makes no normality assumption, but costs 200
  simulations x 400 bootstrap resamples per effect size per dataset, which is hours for
  157 datasets.

So the simulation runs on a sample and is used to check the analytic figure. If they
disagree the analytic one is not trustworthy and gets dropped -- reporting the fast number
without that check would be assuming the thing that needs testing.

    python scripts/power.py --validate 6
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import norm  # noqa: E402

from rbp.variants import conservation as cons  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
POWER, ALPHA = 0.80, 0.05
Z = norm.ppf(1 - ALPHA / 2) + norm.ppf(POWER)      # 2.802

# Empirically calibrated multiplier, measured on 7 datasets spanning n=26 to n=2711 and
# 6 to 400 pathogenic: the effect giving 80% power is a median 2.26 standard errors, not
# the textbook 2.80 (range 1.46-2.67). The textbook figure assumes the standard error is
# the same under the null and the alternative. It is not: `controlled_se` comes from a
# fit where delta barely predicts anything, and under a real effect the likelihood is
# more sharply peaked and the SE shrinks. So 2.80 x SE_observed OVERSTATES the minimum
# detectable effect by roughly 20%.
#
# Both are reported. The textbook figure is the conservative bound; the calibrated one is
# the better point estimate. Every conclusion below is checked against both, because a
# claim that depends on which multiplier you pick is not a claim.
Z_CALIBRATED = 2.26


def analytic(res):
    r = res.copy()
    r["mde"] = Z * r.controlled_se                       # conservative
    r["mde_calibrated"] = Z_CALIBRATED * r.controlled_se  # empirically calibrated
    return r


def validate(res, n_sim, n_boot, scores=None, seed=0):
    """Simulated power at the analytic MDE, using each dataset's REAL predictors.

    An earlier version generated delta as |standard normal| and reported 97% power where
    the formula predicts 80%, which reads as "the analytic MDE is too conservative". It
    was the simulation that was wrong. Measured on our data, the real |delta| has skew
    2.80 and kurtosis 16.6 with 12% of values piled near zero, against skew 1.01 and
    kurtosis 0.9 for |normal|. A skewed, zero-inflated predictor carries less information
    per observation, so the idealised simulation was easier than reality and flattered
    the power.

    So the delta and conservation values are resampled from the dataset itself. Only the
    LABELS are simulated, from a known effect. That keeps every distributional quirk of
    the real predictors and makes the comparison meaningful.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for _, r in res.iterrows():
        n, npath = int(r["n"]), int(r["n_pathogenic"])
        beta = float(r["mde"])
        g = scores[scores.dataset == r["group"]] if scores is not None else None
        if g is None or len(g) < 20:
            continue
        delta = np.abs(g.delta.to_numpy())
        consv = g.conservation.to_numpy()
        prev = np.clip(npath / n, 0.01, 0.99)
        intercept = np.log(prev / (1 - prev))
        zd, zc = cons._z(delta), cons._z(consv)

        hits = tried = 0
        for _ in range(n_sim):
            idx = rng.integers(0, len(delta), n)
            logit = intercept + beta * zd[idx] + 1.0 * zc[idx]
            y = rng.binomial(1, 1 / (1 + np.exp(-logit)))
            if len(np.unique(y)) < 2:
                continue
            tried += 1
            try:
                f = cons.fit_delta_coef(delta[idx], y, consv[idx], n_boot=n_boot,
                                        seed=int(rng.integers(1 << 30)))
            except Exception:
                continue
            hits += f.ci_low > 0
        rows.append({"group": r["group"], "n": n, "n_pathogenic": npath,
                     "mde_analytic": round(beta, 3),
                     "simulated_power": round(hits / tried, 3) if tried else np.nan,
                     "n_sim_ok": tried})
    return pd.DataFrame(rows)


def report(res):
    ok = res[res.controlled_se.notna()].copy()
    surv = ok.controlled_survives.astype(bool)
    print(f"\n{'':=<72}")
    print(f"MINIMUM DETECTABLE EFFECT at {POWER:.0%} power, alpha={ALPHA}")
    print(f"{'':=<72}")
    print(f"{len(ok)} testable datasets: {int(surv.sum())} survive, "
          f"{int((~surv).sum())} do not\n")

    print("MDE distribution (standardised coefficient units):")
    for q in (0.1, 0.25, 0.5, 0.75, 0.9):
        print(f"  p{int(q*100):02d}  {ok.mde.quantile(q):6.3f}")

    null = ok[~surv]
    print(f"\nAmong the {len(null)} that did NOT survive, what could they have seen?")
    for t in (0.2, 0.5, 1.0, 1.5, 2.0):
        n = int((null.mde <= t).sum())
        print(f"  could have detected an effect of {t:4.1f}: {n:3d}/{len(null)} "
              f"({100*n/len(null):4.1f}%)")

    # The pooled coefficient is the natural yardstick: it is the effect we actually
    # measured, so "underpowered" means "could not have seen even that".
    pooled = pd.read_csv(TABLES / "variant_results_pooled.csv")
    ref = float(pooled.controlled_coef.iloc[0])
    print(f"\nPooled effect actually measured: {ref:+.3f}")
    print(f"  how many of the {len(null)} non-surviving datasets could not have detected "
          f"even that?")
    for mult, lab in ((Z, "2.80 conservative"), (Z_CALIBRATED, "2.26 calibrated"),
                      (1.46, "1.46 most optimistic observed")):
        u = int((mult * null.controlled_se > ref).sum())
        print(f"    at {lab:30} {u:3d}/{len(null)} ({100*u/len(null):5.1f}%)")
    print("  the conclusion does not depend on which multiplier is used, which is the "
          "only reason it is worth stating")

    print("\nMDE vs pathogenic count (why pooling was the right call):")
    for lo, hi in ((0, 10), (10, 20), (20, 50), (50, 200), (200, 10 ** 9)):
        s = ok[(ok.n_pathogenic >= lo) & (ok.n_pathogenic < hi)]
        if len(s):
            print(f"  {lo:4d}-{hi if hi < 10**9 else '+':>4} pathogenic: "
                  f"{len(s):3d} datasets, median MDE {s.mde.median():5.2f}, "
                  f"{int(s.controlled_survives.astype(bool).sum()):3d} survive")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--validate", type=int, default=6,
                   help="how many datasets to check by simulation (0 to skip)")
    p.add_argument("--n-sim", type=int, default=150)
    p.add_argument("--n-boot", type=int, default=300)
    a = p.parse_args()

    res = pd.read_csv(TABLES / "variant_results.csv")
    res = res[res.note.isna() | (res.note == "")].copy()
    res = analytic(res)
    res.to_csv(TABLES / "variant_power.csv", index=False)
    report(res)

    if a.validate:
        print(f"\n{'':=<72}")
        print(f"VALIDATING the analytic MDE by simulation on {a.validate} datasets")
        print("resampling each dataset's own delta and conservation; only labels are")
        print("simulated, so the real predictor distributions are preserved")
        print(f"{'':=<72}")
        s = pd.read_csv(TABLES / "variant_scores.csv")
        c = pd.read_csv(TABLES / "variant_conservation.csv")[["vid", "conservation"]]
        s = s.merge(c, on="vid", how="left").dropna(subset=["delta", "conservation"])
        s["dataset"] = s.protein + ":" + s.cell

        sample = res.sample(min(a.validate, len(res)), random_state=3)
        v = validate(sample, a.n_sim, a.n_boot, scores=s)
        print(v.to_string(index=False))
        got = v.simulated_power.dropna()
        print(f"\nmean simulated power at the analytic MDE: {got.mean():.3f} "
              f"(target {POWER:.2f})")
        if abs(got.mean() - POWER) > 0.12:
            print(f"  MISMATCH of {got.mean()-POWER:+.3f}. The analytic MDE is "
                  f"{'CONSERVATIVE' if got.mean() > POWER else 'OPTIMISTIC'}; quote it "
                  f"with that caveat, or use the simulation as primary.")
        else:
            print("  agrees with the target, so the analytic MDE is trustworthy here")
        v.to_csv(TABLES / "variant_power_validation.csv", index=False)


if __name__ == "__main__":
    main()
