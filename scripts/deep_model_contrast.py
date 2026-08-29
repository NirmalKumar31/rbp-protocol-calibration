"""R1g: does the protocol contrast survive a model that is not a bag of k-mers?

    python scripts/deep_model_contrast.py --store ../rbp-store
    python scripts/deep_model_contrast.py --from-cache        # summary from the per-dataset table

THE CLAIM UNDER TEST. R1 measures the NESTED contribution of a 4-mer logistic over a
19-feature composition baseline, and finds it more than twice as large under dinucleotide
matching as under GC matching. The paper's sharpest limitation is that every number in it
comes from one model class. This script runs the identical decomposition on a 3-layer CNN
and on a 19.7M-parameter fine-tuned SpliceBERT.

WHAT IS COMPARED, AND WHY IT IS FAIR. Both arms use the same folds, the same seed, the same
hyperparameters, the same code path and full datasets. The only difference is how the
negative windows were chosen. Nothing here is capped, subsampled or early-stopped
differently between arms, because any of those would confound protocol with training.

THE COMPRESSION CORRECTION IS NOT OPTIONAL. A nested AUROC gain is bounded above by
1 - baseline, and the GC arm's composition baseline is much the higher of the two (0.783
against 0.627), so part of any contrast is arithmetic rather than protocol. The transplant
family from scale_check.py is reproduced here verbatim -- both directions, both links --
because reporting the favourable member alone is question-begging.

THE COMPOSITION BASELINE IS SHARED. Within one arm the reduced model does not depend on
which model score is being added, so it is fitted once per (dataset, arm) and reused. That
is also a check: computing it per model gave bitwise identical values, and this script
asserts it rather than trusting it.
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.eval.nested import gain_over_composition  # noqa: E402
from rbp.utils import panel as panelmod  # noqa: E402

TABLES = ROOT / "results" / "tables"
EVIDENCE = ROOT / "data" / "evidence" / "scores"
MODELS = ["cnn", "splicebert"]
# A dataset whose scored rows fall below this share of its windows is dropped rather than
# analysed on a silent subset. 0.99 keeps every real case (worst observed 0.9986) and would
# still catch a fold that failed to upload, which costs 20%.
MIN_COVERAGE = 0.99
R2 = np.sqrt(2.0)


def log(m):
    print(m, flush=True)


def dprime(a):
    return R2 * norm.ppf(np.clip(a, 1e-6, 1.0 - 1e-6))


def auroc(d):
    return norm.cdf(d / R2)


def logit(a):
    a = np.clip(a, 1e-6, 1.0 - 1e-6)
    return np.log(a / (1.0 - a))


def expit(d):
    return 1.0 / (1.0 + np.exp(-d))


def oof(root, cell, protein, model, folds=5):
    """One out-of-fold score per row, pooled over the folds, or None if incomplete.

    A dataset missing even one fold is dropped rather than pooled from four, because a
    partial pool silently changes which rows the AUROC is computed on.
    """
    parts = []
    for f in range(folds):
        p = root / cell / protein / model / f"fold{f}" / "scores.tsv.gz"
        if not p.exists():
            return None
        parts.append(pd.read_csv(p, sep="\t"))
    s = pd.concat(parts, ignore_index=True)
    if s.id.duplicated().any():
        raise ValueError(f"{protein}:{cell} {model}: a row is scored twice; "
                         "the folds do not partition the data")
    return s[["id", "score"]]


def arm_roots(store):
    """Where each arm's windows and model scores live.

    The dinucleotide arm was scored on Modal months ago and its per-window scores are
    committed to the repository. The GC arm is produced by cloud/modal/modal_gc_sweep.py
    into the local store. Different provenance, identical format.
    """
    store = Path(store)
    return {
        "gc": (store / "processed" / "gc", store / "runs" / "gc"),
        "dn": (store / "processed" / "dinuc", EVIDENCE),
    }


def per_dataset(store, datasets):
    roots = arm_roots(store)
    out = []
    for i, ds in enumerate(datasets, 1):
        protein, cell = ds.split(":")
        row = {"dataset": ds, "protein": protein, "cell": cell}
        ok = True
        for arm, (dataroot, scoreroot) in roots.items():
            f = dataroot / cell / protein / "dataset.tsv"
            if not f.exists():
                ok = False
                break
            d = pd.read_csv(f, sep="\t")
            base = None
            for model in MODELS:
                s = oof(scoreroot, cell, protein, model)
                if s is None:
                    ok = False
                    break
                m = d.merge(s, on="id", how="inner")
                # COVERAGE IS MEASURED, NOT ASSUMED, and it is not always 1.
                #
                # The dinucleotide arm was scored on Modal and its per-window scores are
                # committed here, but the windows themselves were regenerated afterwards.
                # Negative matching is a stochastic search, so 18 of the 94 datasets differ
                # from the scored set by a handful of rows in each direction -- 172 rows in
                # 307,430, or 0.06%. The GC arm cannot have this problem: it is scored from
                # the very dataset.tsv in this store.
                #
                # An inner join is the right handling, because composition features need the
                # sequence and a score with no window is unusable either way. What is not
                # right is doing it silently, so the fraction is recorded per dataset and
                # gated in the summary.
                cov = len(m) / len(d)
                if cov < MIN_COVERAGE:
                    log(f"  SKIP {ds} {model} {arm}: only {len(m)} of {len(d)} rows scored")
                    ok = False
                    break
                row[f"coverage_{arm}"] = min(row.get(f"coverage_{arm}", 1.0), cov)
                g = gain_over_composition(m.seq_rna.values, m.score.values,
                                          m.label.values, m.fold.values)
                # The reduced model is the same fit whichever score is being added to it.
                if base is None:
                    base = g.auroc_composition
                elif abs(base - g.auroc_composition) > 1e-12:
                    raise ValueError(f"{ds} {arm}: composition baseline moved between "
                                     f"models ({base} vs {g.auroc_composition})")
                row[f"comp_{arm}"] = g.auroc_composition
                row[f"{model}_full_{arm}"] = g.auroc_with_score
                row[f"{model}_gain_{arm}"] = g.delta
                row[f"{model}_se_{arm}"] = (g.delta_ci_high - g.delta_ci_low) / (2 * 1.959964)
            if not ok:
                break
            row[f"n_{arm}"] = g.n
        if not ok:
            continue
        out.append(row)
        log(f"[{i:3d}/{len(datasets)}] {ds:18s} "
            + "  ".join(f"{m} {row[f'{m}_gain_dn'] - row[f'{m}_gain_gc']:+.4f}"
                        for m in MODELS))
    return pd.DataFrame(out)


def transplant(m, model):
    """The compression/protocol split, exactly as scale_check.py computes it for the k-mer."""
    q = {}
    for a in ("gc", "dn"):
        m[f"dcomp_{a}"] = dprime(m[f"comp_{a}"])
        m[f"dfull_{a}"] = dprime(m[f"{model}_full_{a}"])
        m[f"dd_{a}"] = m[f"dfull_{a}"] - m[f"dcomp_{a}"]
        m[f"lc_{a}"] = logit(m[f"comp_{a}"])
        m[f"ld_{a}"] = logit(m[f"{model}_full_{a}"]) - m[f"lc_{a}"]
    gain_gc, gain_dn = m[f"{model}_gain_gc"], m[f"{model}_gain_dn"]
    pred_dn = auroc(m.dcomp_dn + m.dd_gc) - auroc(m.dcomp_dn)
    pred_gc = auroc(m.dcomp_gc + m.dd_dn) - auroc(m.dcomp_gc)
    lpred_dn = expit(m.lc_dn + m.ld_gc) - expit(m.lc_dn)
    lpred_gc = expit(m.lc_gc + m.ld_dn) - expit(m.lc_gc)
    q["nested_gc"] = gain_gc.mean()
    q["nested_dn"] = gain_dn.mean()
    q["contrast_auroc"] = (gain_dn - gain_gc).mean()
    q["contrast_dprime"] = (m.dd_dn - m.dd_gc).mean()
    q["contrast_scale_only"] = (pred_dn - gain_gc).mean()
    q["contrast_protocol"] = (gain_dn - pred_dn).mean()
    q["contrast_protocol_reverse"] = (pred_gc - gain_gc).mean()
    q["contrast_protocol_logit"] = (gain_dn - lpred_dn).mean()
    q["contrast_protocol_logit_reverse"] = (lpred_gc - gain_gc).mean()
    return q


def summarise(d, n_boot=2000, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for model in MODELS:
        m = d.copy()
        obs = transplant(m, model)
        idx = rng.integers(0, len(d), size=(n_boot, len(d)))
        boots = {k: [] for k in obs}
        for b in idx:
            for k, v in transplant(d.iloc[b].reset_index(drop=True), model).items():
                boots[k].append(v)
        for k, v in obs.items():
            lo, hi = np.percentile(boots[k], [2.5, 97.5])
            rows.append({"model": model, "quantity": k, "value": v,
                         "ci_low": lo, "ci_high": hi, "n": len(d)})
        pos = (d[f"{model}_gain_dn"] - d[f"{model}_gain_gc"] > 0).sum()
        rows.append({"model": model, "quantity": "contrast_positive_datasets",
                     "value": int(pos), "ci_low": np.nan, "ci_high": np.nan, "n": len(d)})
        # The protocol effect is a RANGE across the transplant family, never one member.
        fam = [obs[k] for k in ("contrast_protocol", "contrast_protocol_reverse",
                                "contrast_protocol_logit", "contrast_protocol_logit_reverse")]
        rows.append({"model": model, "quantity": "protocol_effect_min", "value": min(fam),
                     "ci_low": np.nan, "ci_high": np.nan, "n": len(d)})
        rows.append({"model": model, "quantity": "protocol_effect_max", "value": max(fam),
                     "ci_low": np.nan, "ci_high": np.nan, "n": len(d)})
    for arm in ("gc", "dn"):
        col = f"coverage_{arm}"
        rows.append({"model": "-", "quantity": f"min_row_coverage_{arm}",
                     "value": float(d[col].min()) if col in d else 1.0,
                     "ci_low": np.nan, "ci_high": np.nan, "n": len(d)})
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--from-cache", action="store_true",
                   help="rebuild the summary from the committed per-dataset table")
    p.add_argument("--n-boot", type=int, default=2000)
    p.add_argument("--models", default=None,
                   help="comma-separated subset, for a partial sweep")
    a = p.parse_args()
    warnings.filterwarnings("ignore")
    if a.models:
        global MODELS
        MODELS[:] = [m.strip() for m in a.models.split(",")]

    per = TABLES / "deep_contrast_per_dataset.csv"
    if a.from_cache:
        d = pd.read_csv(per)
    else:
        paired = set(pd.read_csv(TABLES / "rehearsal_binding_gc.csv").dataset)
        paired &= set(pd.read_csv(TABLES / "rehearsal_binding_dinuc.csv").dataset)
        d = per_dataset(a.store, sorted(paired))
        if d.empty:
            sys.exit("no dataset has both arms scored for both models yet")
        d.to_csv(per, index=False)

    log(f"\n=== R1g: deep-model contrast, n = {len(d)} of 94 paired datasets ===")
    s = summarise(d, n_boot=a.n_boot)
    s.to_csv(TABLES / "deep_contrast.csv", index=False)
    for model in MODELS:
        g = s[s.model == model].set_index("quantity")
        log(f"\n{model}")
        log(f"  nested contribution   gc {g.loc['nested_gc','value']:+.4f}   "
            f"dinuc {g.loc['nested_dn','value']:+.4f}")
        log(f"  contrast              {g.loc['contrast_auroc','value']:+.4f} "
            f"[{g.loc['contrast_auroc','ci_low']:+.4f}, "
            f"{g.loc['contrast_auroc','ci_high']:+.4f}]  "
            f"positive in {int(g.loc['contrast_positive_datasets','value'])}/{len(d)}")
        log(f"  protocol effect       {g.loc['protocol_effect_min','value']:+.4f} to "
            f"{g.loc['protocol_effect_max','value']:+.4f}  (across two directions, two links)")
    log(f"\nwrote {per.name} and deep_contrast.csv")


if __name__ == "__main__":
    main()
