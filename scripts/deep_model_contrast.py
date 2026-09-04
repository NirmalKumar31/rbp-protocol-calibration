"""R1g: does the protocol contrast survive a model that is not a bag of k-mers?

    python scripts/deep_model_contrast.py --store ../rbp-store
    python scripts/deep_model_contrast.py --from-cache        # summary from the per-dataset table

THE CLAIM UNDER TEST. R1 measures the NESTED contribution of a 4-mer logistic over a
19-feature composition baseline, and finds it more than twice as large under dinucleotide
matching as under GC matching. The paper's sharpest limitation is that every number in it
comes from one model class. This script runs the identical decomposition on a 3-layer CNN
and on a 19.7M-parameter fine-tuned SpliceBERT.

WHAT IS COMPARED, AND WHY IT IS FAIR. Both arms use the same seed, the same hyperparameters,
the same code path and full datasets, and the only intended difference is how the negative
windows were chosen. Nothing here is capped, subsampled or early-stopped differently between
arms, because any of those would confound protocol with training.

FOLD PROVENANCE, and it is not uniform across arms. For **20 of the 94 dinucleotide-arm datasets** the
committed CNN and SpliceBERT scores were produced under a stratified random partition rather
than config/folds.tsv's chromosome grouping -- fold SIZES preserved, so invisible to any count
check, but up to 23 chromosomes per fold and up to 44.5% of rows having a same-strand neighbour
within 1 kb in a different fold. The GC arm is clean, 94 of 94, and the k-mer is refit here on
the study folds for every dataset, so it is a clean internal control.

Audited and gated by scripts/fold_integrity.py. Dropping the 20 moves the contrasts by at most
0.0038, inside the protein-clustered half-widths, so this is a disclosed limitation rather than
a correction to the claim -- but the CNN's dinucleotide-arm gain is affected and must not be
quoted without it. **Do not restore the "same folds" sentence without rerunning the sweep.**

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

from sklearn.metrics import roc_auc_score  # noqa: E402

from rbp.eval.baseline import oof_scores as kmer_oof  # noqa: E402
from rbp.eval.nested import gain_over_composition  # noqa: E402

TABLES = ROOT / "results" / "tables"
# BOTH ARMS' PER-WINDOW SCORES ARE COMMITTED, and that is what makes this table checkable.
# Until they were, deep_contrast_per_dataset.csv was a terminal artifact: hand-committed, not
# regenerable from anything in the repo, and tied by no assertion to the panel, the hardware
# or the 940 score files. A referee forged it end to end and passed 314/314. 16 MB closes it.
EVIDENCE = ROOT / "data" / "evidence" / "scores"
EVIDENCE_GC = ROOT / "data" / "evidence" / "scores_gc"
EVIDENCE_NEG2 = ROOT / "data" / "evidence" / "scores_neg2"
# The ladder, cheapest first. kmer is refitted in-script on the same rows; the other two are
# read from the per-window scores their sweeps wrote.
MODELS = ["kmer", "cnn", "splicebert"]
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
    roots = {
        "gc": (store / "processed" / "gc",
               EVIDENCE_GC if EVIDENCE_GC.exists() else store / "runs" / "gc"),
        "dn": (store / "processed" / "dinuc", EVIDENCE),
    }
    # The bias-aware arm. Unlike the dinucleotide arm, all 188 of these score sets were
    # verified chromosome-grouped before use (scripts/fold_integrity.py).
    neg2_scores = EVIDENCE_NEG2 if EVIDENCE_NEG2.exists() else store / "runs" / "neg2"
    if (store / "processed" / "neg2").exists() and neg2_scores.exists():
        roots["neg2"] = (store / "processed" / "neg2", neg2_scores)
    return roots


def per_dataset(store, datasets):
    """One row per dataset: every model's nested contribution in both arms, same rows.

    THE COMMON ROW SET IS TAKEN FIRST, BEFORE ANYTHING IS FITTED. Within one arm the three
    models do not all cover the same windows: the k-mer is refitted here so it covers every
    row, while the CNN and SpliceBERT scores were written against the window set as it stood
    when their sweep ran, and 18 of the 94 datasets have since drifted by a handful of rows.

    Fitting each model on whatever it happens to cover would make the composition baseline
    move between models -- observed at AQR:HepG2, 0.59675077 against 0.59675949 -- and a
    ladder whose rungs are measured on different rows is not a ladder. Intersecting first
    costs 0.06% of the rows and makes "same rows, same folds, same estimator" literally true.
    """
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

            # 1. every model's ids, before any fitting
            ids = set(d.id)
            got = {}
            for model in MODELS:
                if model == "kmer":
                    continue            # refitted below, on the common set
                s = oof(scoreroot, cell, protein, model)
                if s is None:
                    ok = False
                    break
                got[model] = s
                ids &= set(s.id)
            if not ok:
                break

            cov = len(ids) / len(d)
            if cov < MIN_COVERAGE:
                log(f"  SKIP {ds} {arm}: only {len(ids)} of {len(d)} rows common to all models")
                ok = False
                break
            row[f"coverage_{arm}"] = cov
            dd = d[d.id.isin(ids)].reset_index(drop=True)

            # 2. the k-mer, refitted on exactly those rows and folds
            if "kmer" in MODELS:
                sc, _, _ = kmer_oof(dd.seq_rna.values, dd.label.values, dd.fold.values, k=4)
                got["kmer"] = pd.DataFrame({"id": dd.id.values, "score": sc})

            # 2b. THE ANCHOR COLUMNS. Raw pooled AUROC and row count computed over every
            # committed score row for this (dataset, model, arm) -- deliberately NOT over the
            # intersection, so verify.py can recompute them from data/evidence alone with no
            # window table and no sequences. This is what ties each row of this table to
            # per-window evidence rather than to itself.
            for model in MODELS:
                if model == "kmer":
                    continue
                # The score file's OWN label column, over ALL its rows. Deliberately not
                # joined to the window table and not restricted to the intersection: the
                # point of this column is that verify.py can reproduce it from
                # data/evidence alone, with no window table and no sequences.
                parts = [pd.read_csv(scoreroot / cell / protein / model / f"fold{k}"
                                     / "scores.tsv.gz", sep="\t") for k in range(5)]
                allsc = pd.concat(parts, ignore_index=True)
                row[f"{model}_raw_{arm}"] = float(
                    roc_auc_score(allsc.label.values, allsc.score.values))
                row[f"{model}_nrows_{arm}"] = int(len(allsc))

            # 3. the nested decomposition, one shared reduced model
            base = None
            for model in MODELS:
                m = dd.merge(got[model], on="id", how="inner")
                if len(m) != len(dd):
                    raise ValueError(f"{ds} {model} {arm}: {len(m)} of {len(dd)} common rows")
                g = gain_over_composition(m.seq_rna.values, m.score.values,
                                          m.label.values, m.fold.values)
                if base is None:
                    base = g.auroc_composition
                elif abs(base - g.auroc_composition) > 1e-9:
                    raise ValueError(f"{ds} {arm}: composition baseline moved between "
                                     f"models ({base} vs {g.auroc_composition})")
                row[f"comp_{arm}"] = g.auroc_composition
                row[f"{model}_full_{arm}"] = g.auroc_with_score
                row[f"{model}_gain_{arm}"] = g.delta
                row[f"{model}_se_{arm}"] = (g.delta_ci_high - g.delta_ci_low) / (2 * 1.959964)
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
    # THE LADDER, AS PAIRED DIFFERENCES. The marginal intervals for the k-mer and the CNN
    # overlap slightly, so "0.0398 < 0.0530" on point estimates is not on its own a claim.
    # The datasets are the same 94 in every rung, so the paired difference is the right
    # statistic and it is what the paper reports.
    # Its own resample, rather than whatever `idx` the last model's loop left behind. The
    # draws are exchangeable so the numbers would not move, but a statistic that silently
    # depends on loop order is the kind of thing that is true until someone reorders MODELS.
    lidx = np.random.default_rng(seed + 1).integers(0, len(d), size=(n_boot, len(d)))
    for a, b in (("cnn", "kmer"), ("splicebert", "cnn"), ("splicebert", "kmer")):
        diff = ((d[f"{a}_gain_dn"] - d[f"{a}_gain_gc"])
                - (d[f"{b}_gain_dn"] - d[f"{b}_gain_gc"]))
        bs = np.array([diff.values[i].mean() for i in lidx])
        lo, hi = np.percentile(bs, [2.5, 97.5])
        rows.append({"model": "ladder", "quantity": f"step_{a}_minus_{b}",
                     "value": diff.mean(), "ci_low": lo, "ci_high": hi, "n": len(d)})
        rows.append({"model": "ladder", "quantity": f"step_{a}_minus_{b}_datasets",
                     "value": int((diff > 0).sum()), "ci_low": np.nan, "ci_high": np.nan,
                     "n": len(d)})

    # THE RATIO SCALE, WHERE THE LADDER REVERSES. This is not optional and it is not a
    # sensitivity check: R1b's own rule is that a result whose sign depends on the scale is
    # not a result unless the reversal has a diagnosis. That rule was applied to the log-odds
    # reversal and then not applied to this paper's own newest headline until a referee
    # pointed it out.
    #
    # The additive contrast grows with capacity. The MULTIPLIER does not: it is ~2.4-3.5x for
    # every model class and is significantly SMALLEST for the largest model, because the
    # AUROC ceiling bites harder at SpliceBERT's higher baseline. Restricted to datasets
    # where every model has a positive gain in both arms, since a log-ratio needs that.
    both = np.ones(len(d), bool)
    for m in MODELS:
        both &= (d[f"{m}_gain_gc"] > 0) & (d[f"{m}_gain_dn"] > 0)
    r = d[both].reset_index(drop=True)
    rows.append({"model": "ratio", "quantity": "datasets_positive_both_arms_all_models",
                 "value": int(both.sum()), "ci_low": np.nan, "ci_high": np.nan, "n": len(d)})
    if both.sum() >= 20:
        ridx = np.random.default_rng(seed + 2).integers(0, len(r), size=(n_boot, len(r)))
        lr = {m: np.log(r[f"{m}_gain_dn"] / r[f"{m}_gain_gc"]) for m in MODELS}
        for m in MODELS:
            bs = np.array([lr[m].values[i].mean() for i in ridx])
            lo, hi = np.percentile(bs, [2.5, 97.5])
            rows.append({"model": "ratio", "quantity": f"log_multiplier_{m}",
                         "value": lr[m].mean(), "ci_low": lo, "ci_high": hi, "n": len(r)})
            rows.append({"model": "ratio", "quantity": f"multiplier_{m}",
                         "value": float(np.exp(lr[m].mean())), "ci_low": float(np.exp(lo)),
                         "ci_high": float(np.exp(hi)), "n": len(r)})
        for m in MODELS:
            # The additive contrast ON THE SAME 77 datasets, so the ratio table's two columns
            # are comparable to each other rather than to the 94-dataset panel above.
            av = (r[f"{m}_gain_dn"] - r[f"{m}_gain_gc"])
            rows.append({"model": "ratio", "quantity": f"additive_contrast_{m}",
                         "value": av.mean(), "ci_low": np.nan, "ci_high": np.nan,
                         "n": len(r)})
        for a, b in (("cnn", "kmer"), ("splicebert", "cnn"), ("splicebert", "kmer")):
            dl = lr[a] - lr[b]
            bs = np.array([dl.values[i].mean() for i in ridx])
            lo, hi = np.percentile(bs, [2.5, 97.5])
            rows.append({"model": "ratio", "quantity": f"logstep_{a}_minus_{b}",
                         "value": dl.mean(), "ci_low": lo, "ci_high": hi, "n": len(r)})
            rows.append({"model": "ratio", "quantity": f"logstep_{a}_minus_{b}_datasets",
                         "value": int((dl > 0).sum()), "ci_low": np.nan, "ci_high": np.nan,
                         "n": len(r)})

    for arm in ("gc", "dn"):
        col = f"coverage_{arm}"
        rows.append({"model": "-", "quantity": f"min_row_coverage_{arm}",
                     # NOT `else 1.0`. Defaulting an absent column to perfect coverage means
                     # DELETING the evidence asserts the strongest possible claim about it --
                     # the same shape as the 2026-08 failure where deleting a file passed.
                     # NaN propagates into a failing gate instead.
                     "value": float(d[col].min()) if col in d else float("nan"),
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
    # `check` makes audit_manuscript.py index every CELL of this table rather than only its
    # column aggregates: its rule is "a `check` column or at most 50 rows", and this table has
    # outgrown 50. Without it, interval bounds quoted in the manuscript read as orphans.
    s.insert(0, "check", s.model + "/" + s.quantity)
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
