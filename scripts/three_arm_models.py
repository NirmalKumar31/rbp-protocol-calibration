"""R1v: the three-protocol span, for all three model classes. The gap every referee named.

    python scripts/three_arm_models.py            # needs the window store
    python scripts/three_arm_models.py --from-cache

THE GAP THIS CLOSES. Until the bias-aware arm was swept for the neural models, everything
three-protocol was 4-mer only (the 5.4-fold span, the transform sweep, the baseline
decomposition, the recommendation test) and everything multi-model was two-arm only. So the
paper's headline was a single-model result and the robustness that defends it belonged to
models that had never seen the third protocol. The two halves of the argument never met in one
cell, and three independent reviewers drew exactly that box.

WHAT IT SHOWS. The span is not a property of the cheap model. It is 5.42x for a 4-mer, 7.63x
for a 7,089-parameter CNN and 3.76x for a fine-tuned 19.78M-parameter SpliceBERT: present for
every model class, and largest for the CNN rather than for the largest model, which is the same
non-monotonicity in capacity that the ratio-scale multiplier already showed.

PROVENANCE, and it matters here more than usual. The dinucleotide arm's committed neural scores
include 20 datasets produced under a partition that is not chromosome-grouped (see
fold_integrity.py). The bias-aware scores swept for this analysis were checked before use: all
188 (dataset, model) sets are chromosome-grouped, maximum 5 chromosomes per fold, zero
misaligned. The neg2 column is therefore clean; the dn column carries the caveat it already had.
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

TABLES = ROOT / "results" / "tables"
ARMS = ("dn", "gc", "neg2")
MODELS = ("kmer", "cnn", "splicebert")
N_BOOT = 4000


def log(m):
    print(m, flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--from-cache", action="store_true",
                   help="use the committed per-dataset table instead of rebuilding it")
    a = p.parse_args()
    warnings.filterwarnings("ignore")

    # THIS TABLE HAD NO PRODUCER, and it is gated. `--store` was accepted and never read: the
    # branch was `if a.from_cache or per.exists()`, so whenever the committed table was
    # present -- always -- it was used, and the flag did nothing. That surfaced after the
    # dinucleotide retrain, when this script printed the pre-retrain contributions while
    # deep_model_contrast.py printed the new ones and neither disagreed with itself.
    #
    # It needs no store. Every column here is a projection of deep_contrast_per_dataset.csv,
    # which IS regenerated from the per-window evidence: checked column by column, the GC and
    # bias-aware columns are bit-identical between the two tables and only the retrained
    # dinucleotide columns differ. So it is rebuilt by projection, and the flag now means what
    # it says.
    per = TABLES / "three_arm_models_per_dataset.csv"
    src = TABLES / "deep_contrast_per_dataset.csv"
    if not a.from_cache and src.exists():
        d = pd.read_csv(src)
        cols = ["dataset", "protein", "cell"]
        for arm in ARMS:
            cols.append(f"comp_{arm}")
            cols += [f"{m}_gain_{arm}" for m in MODELS]
        missing = [c for c in cols if c not in d.columns]
        if missing:
            sys.exit(f"{src} lacks {missing}; run deep_model_contrast.py --store first")
        t = d[cols].copy()
        t.to_csv(per, index=False)
        log(f"rebuilt {per.name} from {src.name} ({len(t)} datasets)")
    elif per.exists():
        t = pd.read_csv(per)
    else:
        sys.exit(f"{per} absent and {src} absent; run deep_model_contrast.py --store first")
    if t.empty:
        sys.exit(f"{per} is empty")

    # Protein-clustered, because 15 of the 79 proteins contribute two datasets each.
    prot = t.protein.to_numpy()
    uniq = np.unique(prot)
    members = [np.flatnonzero(prot == q) for q in uniq]
    rng = np.random.default_rng(0)
    draws = [np.concatenate([members[j] for j in rng.integers(0, len(uniq), len(uniq))])
             for _ in range(N_BOOT)]

    out = []

    def add(check, v, note=""):
        v = np.asarray(v, dtype=float)
        b = np.array([v[i].mean() for i in draws])
        lo, hi = np.percentile(b, [2.5, 97.5])
        out.append({"check": check, "value": float(v.mean()), "ci_low": float(lo),
                    "ci_high": float(hi), "n": len(t), "note": note})
        return float(v.mean())

    log(f"=== R1v: three protocols x three model classes, n = {len(t)} datasets, "
        f"{len(uniq)} proteins ===\n")
    for arm in ARMS:
        add(f"composition AUROC, {arm} arm", t[f"comp_{arm}"])
    log("  composition alone:  " + "  ".join(
        f"{a} {t[f'comp_{a}'].mean():.4f}" for a in ARMS))

    log(f"\n  {'model':12s}" + "".join(f"{a:>12s}" for a in ARMS) + f"{'span':>9s}")
    for m in MODELS:
        vals = []
        for arm in ARMS:
            vals.append(add(f"{m} nested contribution, {arm} arm", t[f"{m}_gain_{arm}"]))
        span = max(vals) / min(vals) if min(vals) > 0 else np.inf
        # The span's own interval, bootstrapped as a ratio of means over the same draws.
        b = []
        for i in draws:
            mm = [t[f"{m}_gain_{arm}"].to_numpy()[i].mean() for arm in ARMS]
            b.append(max(mm) / min(mm) if min(mm) > 0 else np.inf)
        b = np.array(b)
        fin = b[np.isfinite(b)]
        lo, hi = np.percentile(fin, [2.5, 97.5])
        out.append({"check": f"{m} three-protocol span", "value": float(span),
                    "ci_low": float(lo), "ci_high": float(hi), "n": len(t),
                    "note": f"{100 * (~np.isfinite(b)).mean():.1f}% non-finite draws"})
        log(f"  {m:12s}" + "".join(f"{v:+12.4f}" for v in vals)
            + f"{span:8.2f}x  [{lo:.2f}, {hi:.2f}]")

    # THE CLAIM: the span is present for every model class, so it is not an artefact of
    # measuring with a bag of k-mers. Gated as a minimum over models rather than as three
    # separate numbers, because that is the sentence the paper makes.
    spans = {m: float(next(r["value"] for r in out
                           if r["check"] == f"{m} three-protocol span")) for m in MODELS}
    out.append({"check": "minimum three-protocol span over model classes",
                "value": float(min(spans.values())), "n": len(t),
                "note": min(spans, key=spans.get)})
    log(f"\n  minimum span over the three model classes: {min(spans.values()):.2f}x "
        f"({min(spans, key=spans.get)})")

    # And the ordering, which is NOT monotone in capacity: the CNN spans widest.
    order = sorted(spans, key=spans.get)
    out.append({"check": "model class with the widest span", "value": float(max(spans.values())),
                "n": len(t), "note": max(spans, key=spans.get)})
    log("  ordering by span: " + " < ".join(f"{m} {spans[m]:.2f}x" for m in order))
    log("  -> widest for the CNN, not the largest model: the span does not grow with capacity,")
    log("     the same non-monotonicity the ratio-scale multiplier shows.")

    # The bias-aware arm is the lowest for every model, which is the reversal generalised.
    lowest = [min(ARMS, key=lambda arm: t[f"{m}_gain_{arm}"].mean()) for m in MODELS]
    out.append({"check": "model classes whose lowest contribution is the bias-aware arm",
                "value": int(sum(x == "neg2" for x in lowest)), "n": len(MODELS)})
    log(f"\n  the bias-aware arm yields the least for {sum(x == 'neg2' for x in lowest)}"
        f"/{len(MODELS)} model classes, and has the highest composition baseline of the three")

    # THE TWO EFFECTS IN THE SAME UNITS, because the Introduction compares them and a fold
    # change and an AUROC difference are not comparable. Ranges of panel means: protocol range
    # holding the model fixed, model-class range holding the protocol fixed. These are
    # DIFFERENCES of means, which no column aggregate reproduces, so they are emitted here
    # rather than left for audit_manuscript.py to fail on.
    log("\n  the two effects in absolute AUROC, as ranges of panel means:")
    for m in MODELS:
        v = [t[f"{m}_gain_{arm}"].mean() for arm in ARMS]
        out.append({"check": f"protocol range of panel means, within {m}",
                    "value": float(max(v) - min(v)), "n": len(t),
                    "note": f"fold {max(v) / min(v):.2f}"})
        log(f"    protocol range within {m:12s} {max(v) - min(v):.4f}  "
            f"(fold {max(v) / min(v):.2f})")
    for arm in ARMS:
        v = [t[f"{m}_gain_{arm}"].mean() for m in MODELS]
        out.append({"check": f"model-class range of panel means, within {arm} arm",
                    "value": float(max(v) - min(v)), "n": len(t),
                    "note": f"fold {max(v) / min(v):.2f}"})
        log(f"    model-class range within {arm:5s} arm  {max(v) - min(v):.4f}  "
            f"(fold {max(v) / min(v):.2f})")

    pd.DataFrame(out).to_csv(TABLES / "three_arm_models.csv", index=False)
    log("\nwrote three_arm_models.csv")


if __name__ == "__main__":
    main()
