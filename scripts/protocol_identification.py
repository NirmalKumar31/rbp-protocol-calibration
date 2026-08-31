"""R1h: is the "protocol effect" identified? Two attacks a referee made, run properly.

    python scripts/protocol_identification.py

WHAT THIS ANSWERS. R1b splits the contrast into AUROC compression and a residual it calls the
protocol effect, by transplanting one arm's d' increment onto the other arm's baseline. A
statistical referee attacked that on two fronts and both attacks are reproduced here rather
than rebutted in prose.

ATTACK 1: THE STATED RANGE IS AN ARBITRARY TRUNCATION. The paper reports "two directions, two
links" and calls the resulting spread a range. But the transplant is defined for any monotone
link, and the choice moves the estimate further than any interval is wide. Six standard links
are run here, both directions, twelve members in total. VERDICT: the paper's range is too
NARROW, not wrong in sign. The effect stays positive under every member for every model. The
manuscript must widen the range it quotes.

ATTACK 2: THE IDENTIFYING ASSUMPTION IS FALSE, AND THIS ONE LANDS. Transplanting an increment
across baselines requires the increment to be baseline-invariant. It is not. Regressing the
within-arm d' increment on the within-arm d' baseline gives a significantly negative slope for
the k-mer and, far more strongly, for SpliceBERT (-0.34, p = 3.4e-07). Part of any such slope
is mechanical -- the baseline appears on both sides, so estimation noise induces a negative
slope of -Var(noise)/Var(baseline) -- and that part is computed here from Hanley-McNeil
standard errors and subtracted before anything is concluded.

Adding the estimated slope back into the transplant is the honest sensitivity analysis:

    increment carried to the other arm  =  increment + slope * (baseline shift)

VERDICT: the k-mer's protocol effect survives it. SpliceBERT's does not. So the claim "the
protocol effect grows with model capacity" is WITHDRAWN, and what remains is that the k-mer's
protocol effect is positive under every specification tried, while the larger model's is
assumption-dependent and must not be quoted as a point estimate.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import linregress, norm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

TABLES = ROOT / "results" / "tables"
MODELS = ("kmer", "cnn", "splicebert")
R2 = np.sqrt(2.0)


def _c(a):
    return np.clip(a, 1e-6, 1.0 - 1e-6)


# Every monotone map from AUROC to an unbounded-ish scale that a referee could reasonably
# propose. The paper used the first two. Restricting to those two is what made the reported
# range look tighter than the question deserves.
LINKS = {
    "probit": (lambda a: R2 * norm.ppf(_c(a)), lambda z: norm.cdf(z / R2)),
    "logit": (lambda a: np.log(_c(a) / (1 - _c(a))), lambda z: 1 / (1 + np.exp(-z))),
    "arcsine": (lambda a: 2 * np.arcsin(np.sqrt(_c(a))),
                lambda z: np.sin(np.clip(z, 0, np.pi) / 2) ** 2),
    "cloglog": (lambda a: np.log(-np.log(1 - _c(a))), lambda z: 1 - np.exp(-np.exp(z))),
    "logerr": (lambda a: -np.log(1 - _c(a)), lambda z: 1 - np.exp(-np.clip(z, 0, 50))),
    "loga": (lambda a: np.log(_c(a)), lambda z: np.exp(np.clip(z, -50, 0))),
}


def hanley_se(auc, n_pos, n_neg):
    """Hanley-McNeil SE of an AUROC. Used only to size the mechanical regression slope."""
    a = _c(auc)
    q1 = a / (2 - a)
    q2 = 2 * a * a / (1 + a)
    return np.sqrt((a * (1 - a) + (n_pos - 1) * (q1 - a * a)
                    + (n_neg - 1) * (q2 - a * a)) / (n_pos * n_neg))


def baseline_slopes(d):
    """Is the d' increment baseline-invariant? Observed slope, and the part that is artefact.

    The baseline enters the increment with a minus sign, so noise in the baseline estimate
    forces a negative slope even under perfect invariance. That attenuation is
    -Var(noise)/Var(baseline) and is computed here rather than waved at.
    """
    f, _ = LINKS["probit"]
    rows = []
    for m in MODELS:
        for arm in ("gc", "dn"):
            base = f(d[f"comp_{arm}"])
            inc = f(d[f"{m}_full_{arm}"]) - base
            r = linregress(base, inc)
            n = d[f"n_{arm}"].to_numpy()
            se_auc = hanley_se(d[f"comp_{arm}"].to_numpy(), n / 2, n / 2)
            # delta method: sd on the probit scale
            se_base = R2 * se_auc / norm.pdf(norm.ppf(_c(d[f"comp_{arm}"].to_numpy())))
            mech = -np.mean(se_base ** 2) / np.var(base, ddof=1)
            rows.append({"model": m, "arm": arm, "slope": r.slope, "p": r.pvalue,
                         "mechanical_slope": mech,
                         "slope_net_of_mechanism": r.slope - mech})
    return pd.DataFrame(rows)


def protocol(d, m, link, direction, slope=0.0):
    """Protocol effect for one model under one link and direction.

    `slope` is the baseline-invariance correction: the increment is carried across arms with
    a term for how much it would have changed from the baseline shift alone.
    """
    f, inv = LINKS[link]
    cg, cd = f(d.comp_gc), f(d.comp_dn)
    ig = f(d[f"{m}_full_gc"]) - cg
    idn = f(d[f"{m}_full_dn"]) - cd
    if direction == "forward":                 # GC increment -> dinuc baseline
        pred = inv(cd + ig + slope * (cd - cg)) - d.comp_dn
        return d[f"{m}_gain_dn"] - pred
    pred = inv(cg + idn + slope * (cg - cd)) - d.comp_gc
    return pred - d[f"{m}_gain_gc"]


def main():
    d = pd.read_csv(TABLES / "deep_contrast_per_dataset.csv")
    rng = np.random.default_rng(0)
    idx = rng.integers(0, len(d), size=(2000, len(d)))

    def ci(series):
        b = np.array([series.values[i].mean() for i in idx])
        return np.percentile(b, [2.5, 97.5])

    sl = baseline_slopes(d)
    sl.to_csv(TABLES / "protocol_baseline_slopes.csv", index=False)
    print("=== ATTACK 2a: is the d' increment baseline-invariant? ===")
    print(sl.to_string(index=False, float_format=lambda x: f"{x:+.4f}"))

    rows = []
    print("\n=== ATTACK 1: the full link family, 6 links x 2 directions ===")
    for m in MODELS:
        vals = []
        for link in LINKS:
            for direction in ("forward", "reverse"):
                v = protocol(d, m, link, direction).mean()
                vals.append(v)
                rows.append({"check": f"{m}/{link}/{direction}", "model": m, "link": link,
                             "direction": direction, "value": v})
        rows.append({"check": f"{m}/link_family_min", "model": m, "link": "-",
                     "direction": "-", "value": min(vals)})
        rows.append({"check": f"{m}/link_family_max", "model": m, "link": "-",
                     "direction": "-", "value": max(vals)})
        print(f"  {m:11s} range over 12 members: {min(vals):+.4f} to {max(vals):+.4f}"
              f"   all positive: {min(vals) > 0}")

    print("\n=== ATTACK 2b: transplant corrected for the baseline slope ===")
    print(f"  {'model':11s} {'uncorrected (probit fwd)':>26s} {'baseline-ADJUSTED':>26s}")
    for m in MODELS:
        b_gc = sl[(sl.model == m) & (sl.arm == "gc")].slope_net_of_mechanism.iloc[0]
        raw = protocol(d, m, "probit", "forward")
        adj = protocol(d, m, "probit", "forward", slope=b_gc)
        lo, hi = ci(adj)
        rows.append({"check": f"{m}/protocol_baseline_adjusted", "model": m, "link": "probit",
                     "direction": "forward", "value": adj.mean(),
                     "ci_low": lo, "ci_high": hi})
        print(f"  {m:11s} {raw.mean():+26.4f} {adj.mean():+.4f} [{lo:+.4f}, {hi:+.4f}]"
              f"{'   <-- INCLUDES ZERO' if lo <= 0 <= hi else ''}")

    out = pd.DataFrame(rows)
    out.to_csv(TABLES / "protocol_identification.csv", index=False)
    print(f"\nwrote protocol_identification.csv and protocol_baseline_slopes.csv")


if __name__ == "__main__":
    main()
