"""R1w: the transform sweep, the baseline decomposition and the recommendation, per model class.

    python scripts/models_by_protocol.py

WHAT THIS FINISHES. Three of the paper's analyses were 4-mer only, and stated as such: the
eight-transform sweep and its exponent family (R1m/R1u), the protocol-versus-baseline
decomposition (R1n), and the test of the recommendation (R1q). They were 4-mer only for one
reason: the neural models had never been run on the bias-aware arm, so there was no
three-protocol row for them. The neg2 sweep supplied that, and these three analyses need
nothing further -- no GPU, no new training, just the committed per-dataset composition and
gain that the sweep already produced.

WHY IT MATTERS. Every one of the three is a claim about the QUANTITY rather than about the
4-mer, and a claim about a quantity that has only been checked with one estimator is a claim
about that estimator. If the 2.00x floor were a property of bag-of-4-mers measurement it would
move for a fine-tuned transformer; if the baseline explained 11% against the protocol's 1% only
for a k-mer, the mechanism would be a k-mer mechanism.

THE TRANSFORMS AND THE ESTIMATORS ARE IMPORTED, NOT REIMPLEMENTED. scale_sweep.TRANSFORMS and
its fold_range are used directly. Reimplementing an estimator "the same way" is what made R1o
measure a different quantity under the paper's name for a month, and the fix there was an
equality assertion; the fix here is not to have a second copy.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from scale_sweep import TRANSFORMS, fold_range  # noqa: E402

TABLES = ROOT / "results" / "tables"
ARMS = ("dn", "gc", "neg2")
MODELS = ("kmer", "cnn", "splicebert")
N_BOOT = 2000


def log(m):
    print(m, flush=True)


def arm_cols(d, model):
    """(composition, full) per arm for one model. full = composition + nested gain."""
    return {a: (d[f"comp_{a}"].to_numpy(),
                d[f"comp_{a}"].to_numpy() + d[f"{model}_gain_{a}"].to_numpy())
            for a in ARMS}


def main():
    argparse.ArgumentParser().parse_args()
    d = pd.read_csv(TABLES / "three_arm_models_per_dataset.csv")
    out = []

    log(f"n = {len(d)} datasets, three protocols x three model classes\n")

    # --- 1. the transform sweep, per model ----------------------------------------------
    log("=== 1. eight transforms: does the floor hold for every model class? ===\n")
    names = list(TRANSFORMS)
    log(f"  {'transform':34s}" + "".join(f"{m:>13s}" for m in MODELS))
    floors = {}
    for name in names:
        fn = TRANSFORMS[name]
        row = []
        for model in MODELS:
            cf = arm_cols(d, model)
            r = fold_range([fn(*cf[a]) for a in ARMS])
            row.append(r)
            out.append({"check": f"fold range, {name}, {model}", "value": float(r),
                        "n": len(d)})
        log(f"  {name:34s}" + "".join(
            f"{v:13.2f}" if np.isfinite(v) else f"{'n/a':>13s}" for v in row))
    for model in MODELS:
        vals = [v["value"] for v in out
                if v["check"].endswith(model) and v["check"].startswith("fold range,")]
        vals = [v for v in vals if np.isfinite(v)]
        floors[model] = min(vals)
        out.append({"check": f"minimum fold range over transforms, {model}",
                    "value": float(min(vals)), "n": len(d)})
    log("\n  floor over the eight: " + "  ".join(f"{m} {floors[m]:.2f}x" for m in MODELS))

    # --- 2. the exponent family, per model ----------------------------------------------
    log("\n=== 2. is a protocol-free exponent available, and is it the same one? ===\n")
    for model in MODELS:
        cf = arm_cols(d, model)

        def rng(p, cf=cf):
            return fold_range([(f - c) / np.maximum(1 - c, 1e-3) ** p for c, f in
                               (cf[a] for a in ARMS)])

        r = minimize_scalar(lambda p: rng(p) if np.isfinite(rng(p)) else 1e9,
                            bounds=(0.2, 6.0), method="bounded")
        out += [{"check": f"equalising exponent, {model}", "value": float(r.x), "n": len(d)},
                {"check": f"fold range at the equalising exponent, {model}",
                 "value": float(r.fun), "n": len(d)}]
        log(f"  {model:12s} argmin p = {r.x:5.3f}  ->  {r.fun:.3f}x")
    ps = [v["value"] for v in out if v["check"].startswith("equalising exponent")]
    out.append({"check": "spread of the equalising exponent across model classes",
                "value": float(max(ps) / min(ps)), "n": len(MODELS)})
    log(f"\n  the exponent differs {max(ps)/min(ps):.2f}x ACROSS MODEL CLASSES on the same "
        f"benchmark,")
    log("  so it is not even a property of the benchmark alone; it is fitted per measurement.")

    # --- 3. baseline versus protocol label, per model ------------------------------------
    log("\n=== 3. protocol label or composition baseline, per model? ===\n")
    log(f"  {'model':12s}{'protocol|baseline':>19s}{'baseline|protocol':>19s}")
    for model in MODELS:
        long = pd.concat([pd.DataFrame({"arm": a, "comp": d[f"comp_{a}"],
                                        "gain": d[f"{model}_gain_{a}"]}) for a in ARMS],
                         ignore_index=True)
        c, y = long.comp.to_numpy(), long.gain.to_numpy()

        def ss(X, y=y):
            b, *_ = np.linalg.lstsq(X, y, rcond=None)
            r = y - X @ b
            return float(r @ r)

        tss = float((y - y.mean()) @ (y - y.mean()))
        Xb = np.column_stack([np.ones(len(c)), c, c ** 2])
        D = pd.get_dummies(long.arm, drop_first=True).to_numpy(dtype=float)
        Xd = np.column_stack([np.ones(len(c)), D])
        inc_p = (ss(Xb) - ss(np.column_stack([Xb, D]))) / tss
        inc_b = (ss(Xd) - ss(np.column_stack([Xd, c[:, None], (c ** 2)[:, None]]))) / tss
        out += [{"check": f"incremental R2 of the protocol label, {model}",
                 "value": float(inc_p), "n": len(d)},
                {"check": f"incremental R2 of the baseline, {model}", "value": float(inc_b),
                 "n": len(d)}]
        log(f"  {model:12s}{100 * inc_p:18.2f}%{100 * inc_b:18.2f}%")
        # and the within-arm gradient, which is the two-family mechanism
        for a in ARMS:
            r_, p_ = spearmanr(d[f"comp_{a}"], d[f"{model}_gain_{a}"])
            out.append({"check": f"within-arm spearman(baseline, gain), {a}, {model}",
                        "value": float(r_), "n": len(d), "note": f"p={p_:.3f}"})
    log("\n  within-arm gradient (the two-family mechanism), per model:")
    log(f"  {'model':12s}" + "".join(f"{a:>10s}" for a in ARMS))
    for model in MODELS:
        vals = [next(v["value"] for v in out
                     if v["check"] == f"within-arm spearman(baseline, gain), {a}, {model}")
                for a in ARMS]
        log(f"  {model:12s}" + "".join(f"{v:+10.3f}" for v in vals))

    # --- 4. the recommendation, per model ------------------------------------------------
    log("\n=== 4. does the two-number report help, per model? ===\n")
    pairs = (("gc", "dn"), ("gc", "neg2"), ("dn", "neg2"))
    log(f"  {'model':12s}{'rank agreement raw -> headroom':>34s}{'pairs improved':>16s}")
    for model in MODELS:
        raw = {a: d[f"{model}_gain_{a}"].to_numpy() for a in ARMS}
        head = {a: d[f"{model}_gain_{a}"].to_numpy()
                / np.maximum(1 - d[f"comp_{a}"].to_numpy(), 1e-3) for a in ARMS}
        better = 0
        deltas = []
        for x, y in pairs:
            r1 = spearmanr(raw[x], raw[y])[0]
            r2 = spearmanr(head[x], head[y])[0]
            better += int(r2 > r1)
            deltas.append(r2 - r1)
            out += [{"check": f"rank agreement raw, {x} vs {y}, {model}", "value": float(r1)},
                    {"check": f"rank agreement headroom, {x} vs {y}, {model}",
                     "value": float(r2)}]
        out.append({"check": f"protocol pairs where rank agreement improves, {model}",
                    "value": int(better), "n": len(pairs)})
        log(f"  {model:12s}{np.mean(deltas):+34.3f}{better:>13d}/3")
    log("\n  -> mean change in rank agreement under the recommended normalisation, and how")
    log("     many of the three protocol pairs improve, for each model class")

    pd.DataFrame(out).to_csv(TABLES / "models_by_protocol.csv", index=False)
    log("\nwrote models_by_protocol.csv")


if __name__ == "__main__":
    main()
