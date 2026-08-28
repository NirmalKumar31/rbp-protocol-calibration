"""R1: is the protocol effect real, or an artefact of the AUROC scale?

THE OBJECTION. R1's headline is a difference of AUROC differences. The score's nested
contribution over composition is +0.0265 under GC-matched negatives and +0.0662 under
dinucleotide-matched negatives, so the standard protocol appears to hide most of it. But the
two arms do not start from the same place: the composition baseline is 0.7827 in the GC arm
and 0.6280 in the dinucleotide arm. AUROC is bounded at 1, so a fixed increment of real
discriminative signal buys a SMALLER AUROC increment when the baseline is already high. The
compression factor between these two baselines is 1.5x. That predicts the observed direction
with no protocol effect whatsoever, which is why the objection has to be answered with a
number rather than a paragraph.

WHY SOMERS' D IS NOT THE ANSWER. D = 2*AUROC - 1, so a nested contribution on the D scale is
exactly twice the same thing on the AUROC scale and the contrast merely doubles. It is a
linear rescaling and cannot diagnose a nonlinearity. Stated here because it is the first fix
that comes to mind and it is worthless.

WHAT ANSWERS IT. Put both arms on a scale that is linear in signal rather than bounded. Under
a binormal model d' = sqrt(2) * Phi^-1(AUROC), and an increment in d' is an increment in
signal regardless of where the baseline sits. Two quantities follow:

  d' contrast     the same comparison on the unbounded scale. If compression were the whole
                  story this would be zero.

  decomposition   transplant the GC arm's d' increment onto the dinucleotide arm's baseline.
                  That predicts what the dinucleotide arm would show if the protocol changed
                  the baseline and nothing else. The gap between that prediction and what the
                  dinucleotide arm actually shows is the protocol effect with compression
                  removed, and it is the honest headline.

THE THIRD SCALE, AND WHY IT DISAGREES. The committed tables also carry `coef`, the Firth
coefficient of the standardised score in the nested fit. On that scale the contrast REVERSES:
+1.063 in the GC arm against +0.686 in the dinucleotide arm. Reported here rather than buried,
because a result whose sign depends on the scale is not a result -- unless the reversal itself
has a diagnosis, and it does. A logistic coefficient is identified only against the latent
residual scale, so coefficients from two fits with different total signal are not comparable
(Mood 2010, Eur Sociol Rev; the same non-collapsibility that voided this project's earlier
conservation analysis). The GC task carries 1.79x the total signal. The fingerprint is
measured below: the between-arm coefficient gap tracks the TOTAL-signal gap at rho +0.52,
and the incremental-value gap at rho +0.07, p = 0.53. It is measuring the wrong thing, and
dividing it by each fit's own total signal restores the sign.

Bootstrap is a single paired resample of DATASETS -- one stream for every quantity, so the
contrasts stay paired. Differencing two independent streams was a real bug in this project
and inflated an interval 2.6x.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
TABLES = ROOT / "results" / "tables"
N_BOOT = 2000
SEED = 0
R2 = np.sqrt(2.0)


def log(m):
    print(m, flush=True)


def dprime(a):
    return R2 * norm.ppf(np.clip(a, 1e-6, 1.0 - 1e-6))


def auroc(d):
    return norm.cdf(d / R2)


def load():
    gc = pd.read_csv(TABLES / "rehearsal_binding_gc.csv")
    dn = pd.read_csv(TABLES / "rehearsal_binding_dinuc.csv")
    m = gc.merge(dn, on="dataset", suffixes=("_gc", "_dn"))
    for a in ("gc", "dn"):
        m[f"dcomp_{a}"] = dprime(m[f"composition_auroc_{a}"])
        m[f"dfull_{a}"] = dprime(m[f"with_score_auroc_{a}"])
        m[f"dd_{a}"] = m[f"dfull_{a}"] - m[f"dcomp_{a}"]
    # what the dinucleotide arm would show if the GC arm's own added signal were moved onto
    # its baseline: protocol changes the starting point, nothing else
    m["pred_dn"] = auroc(m.dcomp_dn + m.dd_gc) - auroc(m.dcomp_dn)
    return m


def quantities(m):
    """Every contrast, as dinucleotide minus GC. Keys are stable; the bootstrap reuses them."""
    return {
        "contrast_auroc": (m.delta_auroc_dn - m.delta_auroc_gc).mean(),
        "contrast_dprime": (m.dd_dn - m.dd_gc).mean(),
        "contrast_logodds": (m.coef_dn - m.coef_gc).mean(),
        "contrast_scale_only": (m.pred_dn - m.delta_auroc_gc).mean(),
        "contrast_protocol": (m.delta_auroc_dn - m.pred_dn).mean(),
        "contrast_logodds_normalised": (m.coef_dn / m.dfull_dn - m.coef_gc / m.dfull_gc).mean(),
        "nested_gc": m.delta_auroc_gc.mean(),
        "nested_dn": m.delta_auroc_dn.mean(),
    }


def main():
    m = load()
    n = len(m)
    obs = quantities(m)

    rng = np.random.default_rng(SEED)
    boots = {k: [] for k in obs}
    for _ in range(N_BOOT):
        idx = rng.integers(0, n, n)
        for k, v in quantities(m.iloc[idx]).items():
            boots[k].append(v)

    def ci(k):
        return np.percentile(boots[k], [2.5, 97.5])

    rows = []

    def add(check, value, k=None, npt=n, note=""):
        lo, hi = ci(k) if k else (np.nan, np.nan)
        rows.append({"check": check, "value": float(value), "ci_low": lo, "ci_high": hi,
                     "n": npt, "note": note})

    # WHICH MODEL PRODUCED THESE NUMBERS. The manuscript called it a "5-mer" four times,
    # including in the one-sentence claim. Both rehearsal tables record k = 4 on all 189 rows,
    # and `cloud_rehearsal.py:257` defaults --k to KMER_K or 4. The 5 came from
    # `config/params.yaml` `cv: k: 5`, which is the FOLD COUNT of the cross-validation. 150
    # gated assertions could not see this, because every one of them checks a value and none
    # checked what produced it. Emitted here so that it is checkable.
    ks = sorted(set(m.k_gc.unique()) | set(m.k_dn.unique()))
    add("k-mer size, both arms", float(ks[0]) if len(ks) == 1 else -1.0,
        note=f"observed {ks}; -1 means the arms disagree")

    # THE R1 TABLE'S OWN CELLS. Printed in the manuscript, computed on the fly, and stored
    # nowhere until now -- the same state the +0.0397 contrast was in. scripts/audit_manuscript.py
    # found six of them still orphaned after that one was fixed.
    for lab, col in (("composition alone", "composition_auroc"),
                     ("score alone", "auroc"),
                     ("composition + score", "with_score_auroc")):
        for arm, tag in (("gc", "GC"), ("dn", "dinuc")):
            add(f"mean AUROC, {lab}, {tag} arm", float(m[f"{col}_{arm}"].mean()),
                note="cell of the R1 table")

    add("nested contribution, GC arm", obs["nested_gc"], "nested_gc")
    add("nested contribution, dinuc arm", obs["nested_dn"], "nested_dn")
    add("CONTRAST, AUROC scale (published headline)", obs["contrast_auroc"], "contrast_auroc",
        note=f"dinuc larger in {int((m.delta_auroc_dn > m.delta_auroc_gc).sum())}/{n}; "
             "this number appeared in no committed table before now")

    add("compression factor, GC vs dinuc baseline",
        float((norm.pdf(m.dcomp_dn / R2) / norm.pdf(m.dcomp_gc / R2)).mean()),
        note="how much more AUROC a fixed signal increment buys in the dinuc arm")
    add("CONTRAST, d-prime scale (unbounded)", obs["contrast_dprime"], "contrast_dprime",
        note=f"dinuc larger in {int((m.dd_dn > m.dd_gc).sum())}/{n}; zero if compression "
             "were the whole story")

    add("contrast attributable to SCALE alone", obs["contrast_scale_only"],
        "contrast_scale_only", note="GC arm's own d' increment moved onto the dinuc baseline")
    add("CONTRAST, protocol effect net of scale", obs["contrast_protocol"],
        "contrast_protocol",
        note=f"THE HONEST HEADLINE; positive in "
             f"{int((m.delta_auroc_dn > m.pred_dn).sum())}/{n}")
    add("scale share of the published contrast",
        obs["contrast_scale_only"] / obs["contrast_auroc"],
        note="fraction of +0.0397 that is the AUROC ceiling, not the protocol")
    add("fraction of the contribution hidden by GC matching, corrected",
        obs["contrast_protocol"] / obs["nested_dn"],
        note="published claim was two-thirds; uncorrected value is "
             f"{obs['contrast_auroc'] / obs['nested_dn']:.3f}")

    # --- the reversal on the coefficient scale, and its diagnosis ----------------------
    add("CONTRAST, log-odds scale (REVERSES)", obs["contrast_logodds"], "contrast_logodds",
        note=f"dinuc larger in only {int((m.coef_dn > m.coef_gc).sum())}/{n}")
    add("total signal ratio, GC over dinuc",
        float(m.dfull_gc.mean() / m.dfull_dn.mean()),
        note="d'(composition+score); logistic coefficients are not comparable across fits "
             "with different total signal")
    rho_t = spearmanr(m.coef_gc - m.coef_dn, m.dfull_gc - m.dfull_dn)
    rho_i = spearmanr(m.coef_gc - m.coef_dn, m.dd_dn - m.dd_gc)
    add("spearman(coef gap, TOTAL-signal gap)", rho_t.statistic,
        note=f"p={rho_t.pvalue:.3g}; the coefficient gap tracks total signal")
    add("spearman(coef gap, INCREMENTAL-value gap)", rho_i.statistic,
        note=f"p={rho_i.pvalue:.3g}; and not the quantity it is supposed to measure")
    add("CONTRAST, log-odds normalised by total signal",
        obs["contrast_logodds_normalised"], "contrast_logodds_normalised",
        note=f"sign restored; dinuc larger in "
             f"{int((m.coef_dn / m.dfull_dn > m.coef_gc / m.dfull_gc).sum())}/{n}")

    out = pd.DataFrame(rows)
    out.to_csv(TABLES / "scale_check.csv", index=False)
    for _, x in out.iterrows():
        c = f" [{x.ci_low:+.4f}, {x.ci_high:+.4f}]" if pd.notna(x.ci_low) else ""
        log(f"  {x.check:52} {x.value:+.4f}{c}")
    log(f"\n  wrote {TABLES / 'scale_check.csv'}")


if __name__ == "__main__":
    main()
