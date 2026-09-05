"""Is the protocol effect in the training data, the evaluation data, or both?

    python scripts/protocol_transport.py --store ../rbp-store
    python scripts/protocol_transport.py --from-cache

THE CONFOUND. Every arm of this study changes negative construction in training and in
evaluation at the same time. A smaller measured contribution under one protocol could therefore
mean the model learned less from those training negatives, or that those evaluation negatives
are a different discrimination problem, or that the composition baseline moved, or any mixture.
The design as published estimates the effect of changing the whole benchmark protocol, which is
a real and useful estimand, and it cannot separate those.

An external review asked for a train-protocol by evaluation-protocol factorial. This is it, for
the 4-mer, where it costs nothing: fit the model on one arm's windows and score another arm's.

WHY THIS IS SOUND ACROSS ARMS. The chromosome-to-fold map is frozen once for all datasets and
all arms, so fold $i$ of the dinucleotide arm and fold $i$ of the GC arm hold the same
chromosomes. A model fitted on arm A's folds $\\neq i$ has therefore seen none of the
chromosomes in arm B's fold $i$, and scoring across arms leaks nothing that scoring within an
arm would not.

THE VOCABULARY IS THE TRAP. A k-mer count matrix is only meaningful against the vectoriser that
built it; refitting on the evaluation arm would produce a different feature space in which the
training arm's coefficients refer to the wrong columns, silently, since the shapes still match.
rbp.eval.baseline says this in as many words about variant scoring, and it applies identically
here. The training arm's vectoriser transforms the evaluation arm's sequences.

WHAT THE TABLE SEPARATES. With the diagonal being the published within-arm result:

  row effect     holding the evaluation arm fixed and varying the training arm isolates what
                 the training negatives did to the fitted model.
  column effect  holding the training arm fixed and varying the evaluation arm isolates what
                 the evaluation negatives did to the measurement, baseline included.

If the column effect carries nearly all of the protocol dependence, then what the protocol
moves is the measurement rather than the model, which is the paper's thesis stated causally
rather than by association.
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.eval.baseline import fit_fold_models  # noqa: E402
from rbp.eval.nested import gain_over_composition  # noqa: E402
from rbp.utils.log import log  # noqa: E402

TABLES = ROOT / "results" / "tables"
ARMS = {"gc": "gc", "dn": "dinuc", "neg2": "neg2"}
K = 4


def transported_scores(train, evalu, k=K):
    """Score `evalu`'s windows out of fold with models fitted on `train`'s windows.

    Fold discipline is preserved across arms: the model used on evaluation fold i is the one
    that never saw training fold i, and the two arms share a chromosome-to-fold map.
    """
    models, vec = fit_fold_models(list(train.seq_rna), train.label.values,
                                  train.fold.values, k=k)
    X = vec.transform(list(evalu.seq_rna))          # the TRAINING arm's vocabulary
    folds = np.asarray(evalu.fold.values)
    out = np.full(len(evalu), np.nan)
    for f in np.unique(folds):
        m = models.get(int(f))
        if m is None:
            continue
        sel = folds == f
        out[sel] = m.decision_function(X[sel])
    return out


def build(store, limit):
    pub = pd.read_csv(TABLES / "three_arm_per_dataset.csv")
    rows = []
    for n, ds in enumerate(list(pub.dataset)[:limit or None], 1):
        protein, cell = ds.split(":")
        d = {}
        ok = True
        for arm, sub in ARMS.items():
            f = Path(store) / "processed" / sub / cell / protein / "dataset.tsv"
            if not f.exists():
                ok = False
                break
            d[arm] = pd.read_csv(f, sep="\t")
        if not ok:
            continue
        rec = {"dataset": ds, "protein": protein, "cell": cell}
        for tr in ARMS:
            for ev in ARMS:
                sc = transported_scores(d[tr], d[ev])
                keep = np.isfinite(sc)
                if keep.sum() < 200 or d[ev].label.values[keep].std() == 0:
                    ok = False
                    break
                g = gain_over_composition(d[ev].seq_rna.values[keep], sc[keep],
                                          d[ev].label.values[keep], d[ev].fold.values[keep])
                rec[f"gain_{tr}_on_{ev}"] = float(g.delta)
                rec[f"model_{tr}_on_{ev}"] = float(g.auroc_with_score)
            if not ok:
                break
        if not ok:
            continue
        rows.append(rec)
        log(f"[{n:3d}/94] {ds:18s} " + "  ".join(
            f"{tr}->{ev} {rec[f'gain_{tr}_on_{ev}']:+.4f}"
            for tr in ("dn",) for ev in ARMS))
    t = pd.DataFrame(rows)
    if t.empty:
        sys.exit("nothing built; refusing to overwrite the committed table")
    return t


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--n", type=int, default=0)
    p.add_argument("--from-cache", action="store_true")
    a = p.parse_args()
    warnings.filterwarnings("ignore")

    per = TABLES / "protocol_transport_per_dataset.csv"
    t = pd.read_csv(per) if a.from_cache else build(a.store, a.n)
    if not a.from_cache:
        if a.n:
            per = per.with_suffix(".partial.csv")
            log(f"  --n {a.n} given: writing {per.name}, not the committed table")
        t.to_csv(per, index=False)

    rng = np.random.default_rng(0)
    prot = t.protein.to_numpy()
    uniq = np.unique(prot)
    members = [np.flatnonzero(prot == q) for q in uniq]
    draws = [np.concatenate([members[j] for j in rng.integers(0, len(uniq), len(uniq))])
             for _ in range(4000)]
    out = []

    def add(check, v, note=""):
        v = np.asarray(v, dtype=float)
        b = np.array([v[i].mean() for i in draws])
        out.append({"check": check, "value": float(v.mean()),
                    "ci_low": float(np.percentile(b, 2.5)),
                    "ci_high": float(np.percentile(b, 97.5)), "n": len(t), "note": note})

    add("datasets", np.full(len(t), len(t)))
    for tr in ARMS:
        for ev in ARMS:
            add(f"contribution, trained on {tr}, evaluated on {ev}", t[f"gain_{tr}_on_{ev}"],
                "diagonal is the published within-arm value" if tr == ev else "")

    # THE DECOMPOSITION. Both are averages of differences taken with the other factor held, so
    # neither is a main effect from a fitted model; they are the marginal moves the design
    # actually supports.
    cols = list(ARMS)
    tr_eff = np.mean([[t[f"gain_{a1}_on_{ev}"].mean() for a1 in cols] for ev in cols], axis=0)
    ev_eff = np.mean([[t[f"gain_{tr}_on_{a2}"].mean() for a2 in cols] for tr in cols], axis=0)
    out.append({"check": "spread across TRAINING arms, evaluation held fixed",
                "value": float(tr_eff.max() - tr_eff.min()), "ci_low": "", "ci_high": "",
                "n": len(t), "note": "what the training negatives did to the fitted model"})
    out.append({"check": "spread across EVALUATION arms, training held fixed",
                "value": float(ev_eff.max() - ev_eff.min()), "ci_low": "", "ci_high": "",
                "n": len(t), "note": "what the evaluation negatives did to the measurement"})
    tot = (tr_eff.max() - tr_eff.min()) + (ev_eff.max() - ev_eff.min())
    if tot > 0:
        out.append({"check": "share of the protocol effect carried by the evaluation arm",
                    "value": float((ev_eff.max() - ev_eff.min()) / tot), "ci_low": "",
                    "ci_high": "", "n": len(t), "note": ""})

    r = pd.DataFrame(out)
    # The SUMMARY needs the same guard as the per-dataset table: a smoke run must not replace
    # a released artefact. Only the per-dataset path had it, which is how a three-dataset run
    # still overwrote the committed summary.
    summary = TABLES / ("protocol_transport.partial.csv" if a.n else "protocol_transport.csv")
    r.to_csv(summary, index=False)
    log("")
    log("  contribution, rows = training arm, columns = evaluation arm")
    log("           " + "".join(f"{ev:>12s}" for ev in cols))
    for tr in cols:
        log(f"    {tr:6s} " + "".join(f"{t[f'gain_{tr}_on_{ev}'].mean():+12.4f}" for ev in cols))
    log("")
    for _, x in r.tail(3).iterrows():
        log(f"  {x['check']:58s} {x['value']:+.4f}")
    log("\nwrote protocol_transport.csv and protocol_transport_per_dataset.csv")


if __name__ == "__main__":
    main()
