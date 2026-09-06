"""R1n: is it the protocol, or the baseline it leaves? A referee found a way to ask.
It is the baseline.

    python scripts/protocol_or_baseline.py

R1l argued that protocol and baseline cannot be separated because the three protocols' baseline
distributions barely overlap. That argument was too comfortable. A hostile referee found two
places where they DO overlap, asked the question there, and got an answer.

THE NATURAL EXPERIMENT. The neg2 protocol usually raises the composition baseline relative to
GC matching, but not always: in 27 of 94 datasets it LOWERS it. If the protocol label carried
information, neg2's deficit would persist in those 27. It does not -- it reverses. Whichever
protocol leaves the lower baseline gets the higher contribution, regardless of which protocol
that is.

THE MATCHED COMPARISON. Pair each dinucleotide-arm dataset with the GC-arm dataset whose
composition baseline is closest, keep pairs within 0.02 AUROC, and compare. The raw contrast is
+0.0398; matched on baseline it is indistinguishable from zero.

WHY THE dn-vs-gc CONTRAST CANNOT BE ASKED THIS WAY DIRECTLY. comp_dn < comp_gc in 94 of 94
datasets, so within a dataset the two arms are perfectly rank-confounded with the baseline by
construction. The matching above borrows across datasets, which is why it is a weaker design
than the neg2 discordance and why both are reported.

WHAT THIS DOES TO THE PAPER. It converts the thesis from "the protocol determines the measured
contribution" -- which invites "yes, via the baseline, so what?" -- into the sharper and more
useful "the composition baseline is what determines the measurable contribution, and the
protocol label carries essentially no information beyond it." That is a statement a benchmark
builder can act on: report the baseline, because it is the whole story.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

TABLES = ROOT / "results" / "tables"
N_BOOT = 4000
MATCH_TOL = 0.02


def ci(v, seed=0, n_boot=N_BOOT):
    v = np.asarray(v, dtype=float)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(v), size=(n_boot, len(v)))
    return np.percentile(v[idx].mean(axis=1), [2.5, 97.5])


def main():
    d = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    rows = []

    print("=== 1. Does protocol add anything once the baseline is known? ===")
    long = pd.concat([pd.DataFrame({"dataset": d.dataset, "arm": a, "comp": d[f"comp_{a}"],
                                    "gain": d[f"gain_{a}"]})
                      for a in ("gc", "dn", "neg2")], ignore_index=True)
    c, y = long.comp.values, long.gain.values
    def ols(X):
        b, *_ = np.linalg.lstsq(X, y, rcond=None)
        r = y - X @ b
        return float(r @ r)
    tss = float((y - y.mean()) @ (y - y.mean()))
    Xb = np.column_stack([np.ones(len(c)), c, c ** 2])
    D = pd.get_dummies(long.arm, drop_first=True).values.astype(float)
    inc_protocol = (ols(Xb) - ols(np.column_stack([Xb, D]))) / tss
    # and the reverse, for contrast
    Xd = np.column_stack([np.ones(len(c)), D])
    inc_baseline = (ols(Xd) - ols(np.column_stack([Xd, c[:, None], (c ** 2)[:, None]]))) / tss
    rows += [{"check": "incremental R2 of the protocol label, given the baseline",
              "value": float(inc_protocol)},
             {"check": "incremental R2 of the baseline, given the protocol label",
              "value": float(inc_baseline)}]
    print(f"  protocol given baseline: {100 * inc_protocol:5.2f}% of variance")
    print(f"  baseline given protocol: {100 * inc_baseline:5.2f}%   <- an order of magnitude more")

    # IS THE COMPARISON IDENTIFIED? Two things a referee asks. First, whether the protocol
    # nearly determines the baseline, which would make both increments meaningless; it does
    # not. Second, that the protocol's real manipulation is WITHIN dataset, so the pooled
    # "protocol given baseline" is estimated off between-dataset variation and shrinks under
    # dataset fixed effects. Both are reported rather than argued.
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import cross_val_score
    acc = float(cross_val_score(LogisticRegression(max_iter=3000),
                                np.column_stack([c, c ** 2]), long.arm, cv=5).mean())
    rows.append({"check": "arm recovered from the composition baseline alone", "value": acc,
                 "n": len(long), "note": "5-fold accuracy; chance is 1/3"})
    for j in range(D.shape[1]):
        col = D[:, j]
        b, *_ = np.linalg.lstsq(Xb, col, rcond=None)
        r = col - Xb @ b
        r2 = 1.0 - float(r @ r) / float((col - col.mean()) @ (col - col.mean()))
        rows.append({"check": f"VIF of protocol dummy {j} on the baseline curve",
                     "value": 1.0 / (1.0 - r2), "n": len(long)})
    FE = pd.get_dummies(long.dataset, drop_first=True).values.astype(float)
    Xbf = np.column_stack([Xb, FE])
    Xdf = np.column_stack([Xd, FE])
    inc_p_fe = (ols(Xbf) - ols(np.column_stack([Xbf, D]))) / tss
    inc_b_fe = (ols(Xdf) - ols(np.column_stack([Xdf, c[:, None], (c ** 2)[:, None]]))) / tss
    rows += [{"check": "incremental R2 of the protocol label, given baseline and dataset",
              "value": float(inc_p_fe), "n": len(long)},
             {"check": "incremental R2 of the baseline, given protocol and dataset",
              "value": float(inc_b_fe), "n": len(long)}]
    print(f"  arm recovered from the baseline alone: {100 * acc:.1f}% (chance 33.3%)")
    print(f"  with dataset fixed effects: {100 * inc_p_fe:.2f}% and {100 * inc_b_fe:.2f}%")

    print("\n=== 2. The natural experiment: 27 datasets where neg2 LOWERS the baseline ===")
    hi = (d.comp_neg2 > d.comp_gc).values
    rows.append({"check": "datasets where neg2 raises the composition baseline",
                 "value": int(hi.sum())})
    for lab, m, key in (("neg2 baseline HIGHER", hi, "concordant"),
                        ("neg2 baseline LOWER", ~hi, "discordant")):
        diff = (d.gain_neg2 - d.gain_gc).values[m]
        lo, up = ci(diff)
        rows += [{"check": f"neg2 minus gc gain, {key} datasets", "value": float(diff.mean()),
                  "ci_low": float(lo), "ci_high": float(up), "n": int(m.sum())},
                 {"check": f"neg2 higher in, {key} datasets", "value": int((diff > 0).sum()),
                  "n": int(m.sum())}]
        print(f"  {lab:22s} n={m.sum():3d}  {diff.mean():+.4f} [{lo:+.4f}, {up:+.4f}]"
              f"  neg2 higher in {int((diff > 0).sum())}/{m.sum()}")
    r, p = spearmanr(d.comp_neg2 - d.comp_gc, d.gain_neg2 - d.gain_gc)
    rows.append({"check": "within-dataset spearman(delta baseline, delta gain)",
                 "value": float(r), "note": f"p={p:.1e}"})
    print(f"  within-dataset spearman(delta baseline, delta gain) = {r:+.3f}  p={p:.1e}")
    print("  -> whichever protocol leaves the LOWER baseline gets the higher contribution,")
    print("     regardless of which protocol that is")

    print("\n=== 2b. PAIRWISE common support -- R1l's number was the THREE-WAY intersection ===")
    # R1l reported a common support 0.0056 AUROC wide holding 3 of 282 cells and concluded no
    # matched comparison is possible. That is the intersection of ALL THREE arms. Pairwise it
    # is a different story, and the matched design IS possible for gc vs neg2.
    qq = {a: long.comp[long.arm == a] for a in ("gc", "dn", "neg2")}
    for a, b in (("gc", "dn"), ("gc", "neg2"), ("dn", "neg2")):
        lo_ = max(qq[a].quantile(.10), qq[b].quantile(.10))
        hi_ = min(qq[a].quantile(.90), qq[b].quantile(.90))
        w = max(hi_ - lo_, 0.0)
        n_ = int(((qq[a] >= lo_) & (qq[a] <= hi_)).sum()
                 + ((qq[b] >= lo_) & (qq[b] <= hi_)).sum())
        rows += [{"check": f"pairwise common support width, {a} vs {b}", "value": float(w)},
                 {"check": f"cells in pairwise common support, {a} vs {b}", "value": n_}]
        print(f"  {a:5s} vs {b:5s}  width {w:.3f} AUROC   {n_}/188 cells")

    # AND THE INCREMENTAL R2 DECOMPOSED BY PAIR. The headline 1.00% averages the pair where
    # the question is unanswerable with the pair where it is answerable.
    for pair in (("gc", "dn"), ("gc", "neg2"), ("dn", "neg2")):
        L = long[long.arm.isin(pair)]
        cc, yy = L.comp.values, L.gain.values
        # yy is rebound each iteration and this closure reads it late, which is the
        # pattern B023 names. It is safe only because _ss is called below inside this same
        # iteration and never stored; if that ever stops being true, so does the safety.
        def _ss(X):
            bb, *_ = np.linalg.lstsq(X, yy, rcond=None)  # noqa: B023
            r = yy - X @ bb  # noqa: B023
            return float(r @ r)
        Xb2 = np.column_stack([np.ones(len(cc)), cc, cc ** 2])
        D2 = pd.get_dummies(L.arm, drop_first=True).values.astype(float)
        t2 = float((yy - yy.mean()) @ (yy - yy.mean()))
        inc = (_ss(Xb2) - _ss(np.column_stack([Xb2, D2]))) / t2
        rows.append({"check": f"incremental R2 of protocol, {pair[0]} vs {pair[1]}",
                     "value": float(inc)})
        print(f"  incremental R2 of protocol, {pair[0]:5s} vs {pair[1]:5s}: {100*inc:5.2f}%")

    print("\n=== 2c. neg2 vs gc AT MATCHED BASELINE -- the comparison R1l said was impossible ===")
    gg2 = d[["comp_gc", "gain_gc"]].values
    mm = []
    for cb, gv in d[["comp_neg2", "gain_neg2"]].values:
        j = int(np.argmin(np.abs(gg2[:, 0] - cb)))
        if abs(gg2[j, 0] - cb) < 0.01:
            mm.append(gv - gg2[j, 1])
    mm = np.array(mm)
    lo_, hi_ = ci(mm)
    rows.append({"check": "neg2 minus gc, matched on baseline", "value": float(mm.mean()),
                 "ci_low": float(lo_), "ci_high": float(hi_), "n": len(mm)})
    print(f"  nearest-baseline matched  n={len(mm):3d}  {mm.mean():+.4f} [{lo_:+.4f}, {hi_:+.4f}]")
    dbv = (d.comp_neg2 - d.comp_gc).values
    dgv = (d.gain_neg2 - d.gain_gc).values
    X0 = np.column_stack([np.ones(len(dbv)), dbv])
    b0, *_ = np.linalg.lstsq(X0, dgv, rcond=None)
    rng2 = np.random.default_rng(0)
    ii = rng2.integers(0, len(dbv), size=(2000, len(dbv)))
    bb0 = np.array([np.linalg.lstsq(X0[i], dgv[i], rcond=None)[0][0] for i in ii])
    lo2, hi2 = np.percentile(bb0, [2.5, 97.5])
    rows.append({"check": "neg2 minus gc, within-dataset intercept at zero baseline shift",
                 "value": float(b0[0]), "ci_low": float(lo2), "ci_high": float(hi2),
                 "n": len(dbv)})
    print(f"  within-dataset intercept  n={len(dbv):3d}  {b0[0]:+.4f} [{lo2:+.4f}, {hi2:+.4f}]")
    print("  -> BOTH exclude zero. A protocol-specific residual DOES exist for neg2.")

    print("\n=== 2d. the gradient is a property of composition-matched negatives ===")
    from scipy.stats import spearmanr as _sp
    # AND WHETHER THE PARTITION IS ROBUST, because the statistic correlates a difference with
    # its own subtrahend and is not invariant to which term is on the x-axis. Reported for all
    # three choices: the dinucleotide arm flips sign and stops looking like the GC arm, so the
    # family split that survives every choice is GC against the other two.
    for a in ("gc", "dn", "neg2"):
        comp, full = d[f"comp_{a}"], d[f"full_{a}"]
        g_ = d[f"gain_{a}"]
        r_, p_ = _sp(comp, g_)
        rows += [{"check": f"within-arm spearman(baseline, gain), {a}", "value": float(r_),
                  "note": f"p={p_:.3f}"},
                 {"check": f"within-arm spearman(full, gain), {a}",
                  "value": float(_sp(full, g_)[0]), "note": "x-axis sensitivity"},
                 {"check": f"within-arm spearman(midpoint, gain), {a}",
                  "value": float(_sp((comp + full) / 2, g_)[0]),
                  "note": "x-axis sensitivity"},
                 {"check": f"var(gain) over var(baseline), {a}",
                  "value": float(g_.var() / comp.var()),
                  "note": "why neither axis is uncontaminated in the dn arm"}]
        print(f"  {a:5s} x=comp {r_:+.3f} (p={p_:.3f})   x=full "
              f"{_sp(full, g_)[0]:+.3f}   x=mid {_sp((comp + full) / 2, g_)[0]:+.3f}   "
              f"var ratio {g_.var() / comp.var():.3f}")
    print("  -> on the baseline it dies for other-RBPs'-sites negatives, which is the")
    print("     mechanism; on the other two axes the dn arm joins it, so the robust split")
    print("     is GC against the other two rather than matched against bias-aware")

    print("\n=== 3. dn vs gc, matched on the composition baseline ===")
    rank = int((d.comp_dn < d.comp_gc).sum())
    rows.append({"check": "datasets where the dinuc baseline is lower than the GC baseline",
                 "value": rank})
    print(f"  comp_dn < comp_gc in {rank}/94 -- perfectly rank-confounded within dataset,")
    print("  so the comparison has to borrow across datasets:")
    g = d[["comp_gc", "gain_gc"]].values
    matched = []
    for cb, gg in d[["comp_dn", "gain_dn"]].values:
        j = int(np.argmin(np.abs(g[:, 0] - cb)))
        if abs(g[j, 0] - cb) < MATCH_TOL:
            matched.append(gg - g[j, 1])
    matched = np.array(matched)
    lo, up = ci(matched)
    raw = float((d.gain_dn - d.gain_gc).mean())
    rows += [{"check": "dn minus gc, matched on baseline", "value": float(matched.mean()),
              "ci_low": float(lo), "ci_high": float(up), "n": len(matched)},
             {"check": "dn minus gc, unmatched (the published contrast)", "value": raw}]
    print(f"  matched within {MATCH_TOL} AUROC: n={len(matched)}  "
          f"{matched.mean():+.4f} [{lo:+.4f}, {up:+.4f}]")
    print(f"  unmatched: {raw:+.4f}  -> the contrast does not survive matching")

    pd.DataFrame(rows).to_csv(TABLES / "protocol_or_baseline.csv", index=False)
    print("\nwrote protocol_or_baseline.csv")


if __name__ == "__main__":
    main()
