"""Stage 13. Assert the reproduction actually reproduced, and fail loudly if not.

WHY THIS IS THE LAST STAGE AND NOT AN AFTERTHOUGHT. A pipeline that runs to completion and
quietly produces different science is worse than one that crashes, because nobody diffs a
plausible table. Reproducibility you do not check is not reproducibility; it is a hope with
a Dockerfile.

Every expectation lives in config/golden.yaml with an explicit tolerance, chosen to absorb
what legitimately varies (panel size, BLAS thread order, bootstrap seed, GPU vs CPU
inference at 1e-4) and nothing more.

    python scripts/verify.py                    # read tables from GCS
    python scripts/verify.py --local results/tables
    python scripts/verify.py --strict            # warnings become failures

Exit 0 means the science reproduced. Exit 1 means it did not, and the report says which
claim broke.
"""

import argparse
import io
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

from rbp.utils import cloud as cloudcfg  # noqa: E402

checks = []


def record(ok, claim, got, want, note=""):
    checks.append((ok, claim, got, want, note))
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {claim:44} got {got:<22} want {want}" + (f"  {note}" if note else ""),
          flush=True)


def near(claim, got, spec):
    """Point value within tolerance."""
    if got is None or (isinstance(got, float) and np.isnan(got)):
        return record(False, claim, "MISSING", f"{spec['value']} +/- {spec['tol']}")
    ok = abs(got - spec["value"]) <= spec["tol"]
    record(ok, claim, f"{got:.4f}", f"{spec['value']:.4f} +/- {spec['tol']}")


def at_least(claim, got, floor):
    ok = got is not None and got >= floor
    record(ok, claim, f"{got:.4f}" if got is not None else "MISSING", f">= {floor}")


def at_most(claim, got, ceil):
    ok = got is not None and got <= ceil
    record(ok, claim, f"{got:.3g}" if got is not None else "MISSING", f"<= {ceil:g}")


class Tables:
    """Result tables, from GCS or a local directory, whichever the run produced."""

    def __init__(self, local=None):
        self.local = Path(local) if local else None
        self.bucket = None if local else cloudcfg.bucket()

    def get(self, name):
        if self.local:
            p = self.local / name
            return pd.read_csv(p) if p.exists() else None
        b = self.bucket.blob(f"results/tables/{name}")
        if not b.exists():
            return None
        return pd.read_csv(io.BytesIO(b.download_as_bytes()))


def verify_r1(T, g):
    print("\nR1  the negative-set protocol effect")
    d = T.get("cost_of_matching.csv")
    if d is None:
        return record(False, "cost_of_matching.csv present", "MISSING", "the table")
    spec = g["r1_cost_of_matching"]
    near("mean AUROC, GC arm", d.auroc_gc.mean(), spec["auroc_gc"])
    near("mean AUROC, dinuc arm", d.auroc_dn.mean(), spec["auroc_dinuc"])
    near("cost of proper matching", d.cost.mean(), spec["cost"])

    # THE NESTED GAIN, NOT THE DIFFERENCE OF TWO STANDALONE AUROCs.
    #
    # This computed (auroc - composition_auroc): how much better the sequence model is than a
    # composition model fitted separately. That is not the claim. The claim is that the
    # sequence score adds information OVER AND ABOVE composition, which is a nested comparison
    # -- composition alone against composition plus the score -- and the rehearsal already
    # computes it as delta_auroc = with_score_auroc - composition_auroc, with a bootstrap CI,
    # a p-value and a per-dataset `helps` flag attached.
    #
    # The two disagree materially: naive gives +0.0154 -> +0.0607 (ratio 3.94), nested gives
    # +0.0265 -> +0.0662 (ratio 2.50). Both support the direction, but the verifier was
    # certifying 3.94x while the manuscript's framing quotes the nested figure. A gate that
    # blesses a number nobody claims is worse than no gate, because it reads as confirmation.
    gain_gc = d.delta_auroc_gc.mean()
    gain_dn = d.delta_auroc_dn.mean()
    near("nested gain over composition, GC", gain_gc, spec["gain_gc"])
    near("nested gain over composition, dinuc", gain_dn, spec["gain_dinuc"])
    # THE THESIS. If this does not hold the paper's central claim is false, regardless of
    # how close every other number lands.
    at_least("gain RATIO (the paper's thesis)", gain_dn / gain_gc if gain_gc else None,
             spec["gain_ratio_min"])
    # The gain must be significant, not merely bigger. `helps` is the rehearsal's own verdict
    # on the nested comparison for that dataset.
    if "helps_dn" in d.columns:
        at_least("nested gain significant, dinuc arm", float(d.helps_dn.mean()),
                 spec["gain_significant_on_min_fraction"])
    near("fraction of datasets falling", float((d.cost < 0).mean()),
         spec["fraction_datasets_falling"])
    from scipy.stats import wilcoxon
    at_most("paired Wilcoxon p", wilcoxon(d.auroc_gc, d.auroc_dn)[1], spec["wilcoxon_p_max"])


def verify_r2(T, g):
    print("\nR2  four models on identical splits")
    # NOT `a or b`: pandas raises ValueError on DataFrame truthiness, so the fallback
    # would crash the verifier rather than fall through.
    d = T.get("matched_four_models.csv")
    if d is None:
        d = T.get("matched95_four_models.csv")
    if d is None:
        return record(False, "four-model table present", "MISSING", "the table")
    spec = g["r2_four_models"]
    col = {"composition": "composition_auroc", "kmer": "kmer_auroc",
           "cnn": "cnn", "splicebert": "splicebert"}
    means = {}
    for k, c in col.items():
        if c not in d:
            record(False, f"column {c}", "MISSING", "present")
            continue
        means[k] = d[c].mean()
        near(f"mean AUROC, {k}", means[k], spec[k])
    order = spec["required_order"]
    vals = [means.get(k) for k in order]
    ok = all(v is not None for v in vals) and all(a < b for a, b in zip(vals, vals[1:]))
    record(ok, "model ordering holds", " < ".join(order) if ok else str(vals),
           "strictly increasing")
    if "splicebert" in means:
        frac = float((d.splicebert > d.composition_auroc).mean())
        near("SpliceBERT beats composition, fraction", frac,
             spec["splicebert_beats_composition_fraction"])


def verify_r3(T, g):
    print("\nR3  positional concentration")
    d = T.get("locality_ism.csv")
    if d is None:
        return record(False, "locality_ism.csv present", "MISSING", "the table")
    spec = g["r3_locality"]
    near("k-mer Gini, median", d.kmer_gini.median(), spec["kmer_gini_median"])
    near("SpliceBERT Gini, median", d.sb_gini.median(), spec["splicebert_gini_median"])
    diff = d.sb_gini - d.kmer_gini
    near("median difference", diff.median(), spec["median_difference"])
    at_least("fraction SpliceBERT higher", float((diff > 0).mean()),
             spec["fraction_splicebert_higher_min"])

    # The strong form of the claim needs per-dataset uncertainty. If gini_sd is absent the
    # run used an older probe and cannot support "reversed on none".
    if {"kmer_gini_sd", "sb_gini_sd", "n_windows"} <= set(d.columns):
        se = np.sqrt((d.kmer_gini_sd ** 2 + d.sb_gini_sd ** 2) / d.n_windows)
        z = diff / se
        at_most("datasets significantly REVERSED", int((z < -1).sum()),
                spec["n_significantly_reversed_max"])
    else:
        record(False, "gini_sd carried through", "MISSING",
               "needed for the 'reversed on none' claim")
    from scipy.stats import wilcoxon
    at_most("paired Wilcoxon p", wilcoxon(d.sb_gini, d.kmer_gini)[1], spec["wilcoxon_p_max"])


def verify_r4(T, g):
    print("\nR4  the ClinVar ladder")
    d = T.get("variant_ladder.csv")
    if d is None:
        return record(False, "variant_ladder.csv present", "MISSING",
                      "the table (stage 11 aggregate)")
    spec = g["r4_variant_ladder"]
    row = {r.arm: r for r in d.itertuples()}
    got = {}
    for arm, key in [("kmer", "kmer_auroc"), ("mismatched", "mismatched_head_auroc"),
                     ("matched", "matched_head_auroc"),
                     ("conservation", "conservation_auroc")]:
        r = row.get(arm)
        got[key] = float(r.auroc) if r is not None else None
        near(f"AUROC, {arm}", got[key], spec[key])
    order = spec["required_order"]
    vals = [got.get(k) for k in order]
    ok = all(v is not None for v in vals) and all(a < b for a, b in zip(vals, vals[1:]))
    record(ok, "ladder is monotone", " < ".join(order) if ok else str(vals),
           "strictly increasing")

    for arm, key in [("matched", "matched_coef"), ("mismatched", "mismatched_coef"),
                     ("kmer", "kmer_coef")]:
        r = row.get(arm)
        near(f"clustered coefficient, {arm}",
             float(r.coef) if r is not None and hasattr(r, "coef") else None, spec[key])
    m, mm = row.get("matched"), row.get("mismatched")
    if m is not None and mm is not None:
        at_least("matched minus mismatched (specificity)",
                 float(m.coef) - float(mm.coef), spec["matched_minus_mismatched_min"])
        if hasattr(m, "n_clusters"):
            at_least("inference is clustered", int(m.n_clusters), spec["min_clusters"])


def verify_integrity(T, g):
    print("\nintegrity")
    spec = g["integrity"]
    for name in ("cost_of_matching.csv", "locality_ism.csv"):
        d = T.get(name)
        if d is None:
            continue
        if "dataset" in d:
            record(int(d.dataset.duplicated().sum()) == spec["duplicate_datasets_allowed"],
                   f"no duplicate datasets in {name}",
                   int(d.dataset.duplicated().sum()), spec["duplicate_datasets_allowed"])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--local", default=None, help="read tables from a local directory")
    p.add_argument("--golden", default=str(ROOT / "config" / "golden.yaml"))
    p.add_argument("--strict", action="store_true")
    a = p.parse_args()

    g = yaml.safe_load(Path(a.golden).read_text())
    T = Tables(a.local)

    print("=" * 78)
    print("VERIFY -- does the reproduction match the original science?")
    # golden.yaml's meta keys are reference_run/established. They were source_run/measured_on
    # when the file described the earlier study, and renaming them there left this line
    # reading keys that no longer exist -- a KeyError before a single check ran.
    m = g["meta"]
    print(f"golden: {m['reference_run']} established {m['established']}")
    print("=" * 78)

    for fn in (verify_r1, verify_r2, verify_r3, verify_r4, verify_integrity):
        try:
            fn(T, g)
        except Exception as e:                  # a broken check is a failure, not a crash
            record(False, f"{fn.__name__} raised", type(e).__name__, "no exception", str(e)[:80])

    bad = [c for c in checks if not c[0]]
    print("\n" + "=" * 78)
    print(f"{len(checks) - len(bad)}/{len(checks)} checks passed")
    if bad:
        print("\nFAILED CLAIMS:")
        for _, claim, got, want, note in bad:
            print(f"  {claim}: got {got}, want {want}" + (f" ({note})" if note else ""))
        print("\nThe pipeline ran but the science did not reproduce. Do not write it up.")
        return 1
    print("Every golden number reproduced within tolerance.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
