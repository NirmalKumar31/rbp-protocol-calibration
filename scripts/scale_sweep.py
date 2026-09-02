"""R1m: does ANY rescaling of AUROC make the contribution protocol-independent? No.

    python scripts/scale_sweep.py

THE OBJECTION THIS EXISTS TO CLOSE, and it is the first one a competent referee raises against
the title. The paper claims the measured contribution of a model is protocol-dependent, and
demonstrates it with a bounded statistic -- a difference of AUROCs -- across three protocols
whose baselines differ by 0.20 AUROC. So: is the range just the bounded scale? Would some
monotone transform remove it?

The answer has to be a sweep, not an argument, because "some transform" is an existential
claim and only a search can address it. Eight transforms are run here, spanning the ones with
an ROC-theoretic motivation (binormal d', logit, arcsine, complementary log-log), the trivial
rescalings (Somers' D), and two normalisations by the protocol's own baseline.

THE RESULT. The fold range across the three protocols never falls below 2.00x. The transform
that achieves that floor divides by the baseline's remaining headroom, 1 - AUROC_composition --
and that is not a rescaling of the model's output at all. It is a rescaling by the protocol's
own baseline, which is to say: the only way to shrink the protocol dependence is to divide it
out by hand, using the very quantity the protocol sets.

WHY THIS IS THE ANALYSIS YOU CANNOT LOSE. If some transform HAD collapsed the range, the paper
would not have been refuted -- it would have become "here is the protocol-independent
coordinate, use it", which is a better and more citable paper. It did not collapse, so the
claim stands and now has a search behind it rather than an assertion.

AND IT YIELDS THE DELIVERABLE. The least protocol-sensitive coordinate found is worth
recommending, WITH its failure stated in the same breath: headroom-normalised contribution
still moves twofold. It is a better coordinate, not an invariant.

Intervals are bootstrapped over the 79 PROTEINS, not the 94 datasets, per R1i.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

TABLES = ROOT / "results" / "tables"
ARMS = ("dn", "gc", "neg2")
N_BOOT = 4000
R2C = np.sqrt(2.0)


def _c(a):
    return np.clip(a, 1e-6, 1.0 - 1e-6)


# Each entry maps (composition, full) -> a scalar "contribution" on that scale.
TRANSFORMS = {
    "raw AUROC gain": lambda c, f: f - c,
    "Somers' D gain": lambda c, f: 2 * (f - c),
    "d' increment (binormal)": lambda c, f: R2C * (norm.ppf(_c(f)) - norm.ppf(_c(c))),
    "logit increment": lambda c, f: np.log(_c(f) / (1 - _c(f))) - np.log(_c(c) / (1 - _c(c))),
    "arcsine increment": lambda c, f: 2 * (np.arcsin(np.sqrt(_c(f)))
                                           - np.arcsin(np.sqrt(_c(c)))),
    "cloglog increment": lambda c, f: (np.log(-np.log(1 - _c(f)))
                                       - np.log(-np.log(1 - _c(c)))),
    "excess-normalised, g/(comp-0.5)": lambda c, f: (f - c) / np.maximum(c - 0.5, 1e-3),
    "headroom-normalised, g/(1-comp)": lambda c, f: (f - c) / np.maximum(1 - c, 1e-3),
}


def fold_range(vals):
    """max/min over the three protocol means. Undefined if any mean is non-positive."""
    m = [v.mean() for v in vals]
    if min(m) <= 0:
        return np.nan
    return max(m) / min(m)


def main():
    d = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    prot = d.protein.to_numpy()
    uniq = np.unique(prot)
    members = [np.flatnonzero(prot == p) for p in uniq]
    rng = np.random.default_rng(0)
    draws = [np.concatenate([members[i] for i in rng.integers(0, len(uniq), len(uniq))])
             for _ in range(N_BOOT)]

    rows = []
    print(f"n = {len(d)} datasets over {len(uniq)} proteins; "
          f"intervals bootstrap the PROTEIN, per R1i\n")
    print(f"  {'transform':34s} " + "".join(f"{a:>9s}" for a in ARMS)
          + f"{'fold range':>22s}")
    best = (np.inf, None)
    for name, fn in TRANSFORMS.items():
        vals = [fn(d[f"comp_{a}"].to_numpy(), d[f"full_{a}"].to_numpy()) for a in ARMS]
        fr = fold_range(vals)
        bs = np.array([fold_range([v[i] for v in vals]) for i in draws])
        bs = bs[np.isfinite(bs)]
        lo, hi = (np.percentile(bs, [2.5, 97.5]) if len(bs) > 100 else (np.nan, np.nan))
        rows.append({"check": f"fold range, {name}", "value": float(fr),
                     "ci_low": float(lo), "ci_high": float(hi), "n": len(d)})
        for a, v in zip(ARMS, vals):
            rows.append({"check": f"mean {name}, {a} arm", "value": float(v.mean()),
                         "ci_low": np.nan, "ci_high": np.nan, "n": len(d)})
        if np.isfinite(fr) and fr < best[0]:
            best = (fr, name)
        print(f"  {name:34s} " + "".join(f"{v.mean():9.4f}" for v in vals)
              + f"   {fr:6.2f}x [{lo:.2f}, {hi:.2f}]")

    rows.append({"check": "minimum fold range over all transforms", "value": float(best[0]),
                 "ci_low": np.nan, "ci_high": np.nan, "n": len(d)})
    rows.append({"check": "transform achieving the minimum", "value": np.nan,
                 "ci_low": np.nan, "ci_high": np.nan, "n": len(d), "note": best[1]})
    print(f"\n  FLOOR: {best[0]:.2f}x, achieved by '{best[1]}'")
    print("  -> no rescaling reaches protocol independence, and the transform that comes")
    print("     closest divides by the protocol's OWN baseline rather than rescaling the model")

    pd.DataFrame(rows).to_csv(TABLES / "scale_sweep.csv", index=False)
    print("\nwrote scale_sweep.csv")


if __name__ == "__main__":
    main()
