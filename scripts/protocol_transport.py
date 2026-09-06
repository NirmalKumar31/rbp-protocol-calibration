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

If the column effect carries most of the protocol dependence, what the protocol moves is
largely the measurement rather than the model. This is DESCRIPTIVE and not causal, for three
reasons that belong next to the number rather than in a later caveat: the off-diagonal cells
evaluate a model on windows drawn under a protocol it was not fitted for, which is a
distribution shift as well as a protocol change; the contribution is computed with the
two-stage estimator whose information route this paper identifies; and the shares depend on how
datasets are weighted, so both weightings are reported.
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

    # A TWO-WAY DECOMPOSITION, WITH AN INTERACTION TERM AND AN INTERVAL.
    #
    # The first version of this reported ratio-of-ranges: it took the range of the three
    # training-arm marginal means (0.0141), the range of the three evaluation-arm marginal means
    # (0.0438), and called 0.0438/(0.0141+0.0438) = 76% "the share carried by the evaluation
    # arm". An audit pointed out that is not a decomposition of anything. It has no interaction
    # term, no uncertainty, and depends on the range functional and on which three arms were
    # chosen. Two ranges summing to a whole is an arithmetic coincidence, not a partition.
    #
    # This is the ordinary balanced two-way sum-of-squares split, computed PER DATASET so the
    # protein-clustered bootstrap can carry it, and reporting the interaction rather than
    # folding it into one of the main effects.
    cols = list(ARMS)
    M = np.stack([[[t[f"gain_{tr}_on_{ev}"].to_numpy(float) for ev in cols] for tr in cols]])
    M = M[0]                                              # 3 x 3 x n_datasets
    mu = M.mean(axis=(0, 1))
    row = M.mean(axis=1) - mu                             # training main effect, 3 x n
    col = M.mean(axis=0) - mu                             # evaluation main effect, 3 x n
    inter = M - M.mean(axis=1)[:, None, :] - M.mean(axis=0)[None, :, :] + mu
    ss_tr = len(cols) * (row ** 2).sum(axis=0)
    ss_ev = len(cols) * (col ** 2).sum(axis=0)
    ss_in = (inter ** 2).sum(axis=(0, 1))
    tot = ss_tr + ss_ev + ss_in
    good = tot > 0
    W = ("mean over datasets of each dataset's own normalised share; equal weight per dataset, "
         "so a dataset with a tiny protocol effect counts as much as a large one")
    add("share of variance from the TRAINING protocol, per-dataset weighting",
        ss_tr[good] / tot[good], W)
    add("share of variance from the EVALUATION protocol, per-dataset weighting",
        ss_ev[good] / tot[good], W)
    add("share of variance from their INTERACTION, per-dataset weighting",
        ss_in[good] / tot[good], W)

    # THE OTHER WEIGHTING, because the choice is not neutral and naming one without the other
    # invites the reader to assume there is only one answer. This decomposes the 3x3 matrix of
    # PANEL MEANS, which weights datasets by how large their effect is. It gives a larger
    # evaluation share and a smaller interaction; both are legitimate and they answer different
    # questions. Reported so the dependence on weighting is visible rather than buried.
    P = M.mean(axis=2)                                    # 3 x 3 panel means
    mp = P.mean()
    rp, cp = P.mean(axis=1) - mp, P.mean(axis=0) - mp
    ip = P - P.mean(axis=1)[:, None] - P.mean(axis=0)[None, :] + mp
    sp_tr, sp_ev, sp_in = 3 * (rp ** 2).sum(), 3 * (cp ** 2).sum(), (ip ** 2).sum()
    tp = sp_tr + sp_ev + sp_in
    for nm, v in (("TRAINING", sp_tr), ("EVALUATION", sp_ev), ("INTERACTION", sp_in)):
        out.append({"check": f"share of variance from the {nm} protocol, panel-mean weighting"
                             if nm != "INTERACTION" else
                             "share of variance from their INTERACTION, panel-mean weighting",
                    "value": float(v / tp), "ci_low": "", "ci_high": "", "n": len(t),
                    "note": "decomposition of the 3x3 matrix of panel means; no interval, "
                            "because it is a function of means rather than a per-dataset value"})

    # The marginal ranges are kept as the descriptive quantities they are, unnormalised, so a
    # reader can see the raw movement without a ratio being read as a partition.
    tr_eff = np.mean([[t[f"gain_{a1}_on_{ev}"].mean() for a1 in cols] for ev in cols], axis=0)
    ev_eff = np.mean([[t[f"gain_{tr}_on_{a2}"].mean() for a2 in cols] for tr in cols], axis=0)
    out.append({"check": "range of training-arm marginal means",
                "value": float(tr_eff.max() - tr_eff.min()), "ci_low": "", "ci_high": "",
                "n": len(t), "note": "descriptive; not a variance component"})
    out.append({"check": "range of evaluation-arm marginal means",
                "value": float(ev_eff.max() - ev_eff.min()), "ci_low": "", "ci_high": "",
                "n": len(t), "note": "descriptive; not a variance component"})

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
