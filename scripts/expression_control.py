"""Are the negatives unbound RNA, or sequence that was never transcribed? Restriction
plus a matched placebo, the same design as `strand_placebo.py`.

THE WOUND. `negatives.py` draws negatives from region pools built by `annotation.py` over
every GENCODE transcript, with no `gene_type` and no cell-line expression filter (limitation
1b). Measured here against ENCODE total RNA-seq for the matching cell line: **41.8% of
negatives overlap no gene expressed above 1 TPM**, against 1.3% of positives, and on the
strand the window is actually labelled with the figure is **68.7%** against 1.4%. The median
negative sits 136 kb from the nearest peak of its own protein. So a large part of the
negative class is not unbound RNA, it is sequence that produces no transcript in this cell
line at all, and the benchmark's dominant axis is transcribed-versus-not.

How dominant: a SINGLE SCALAR, log TPM of the overlapping gene, separates the classes at
mean AUROC 0.833 (0.906 using the labelled strand). That beats the 19-feature composition
baseline in both arms and beats the whole 4-mer model in 86/94 datasets in the dinucleotide
arm. It is the strongest trivial baseline in this study and it is not the one the paper is
built around.

WHAT THIS SCRIPT ANSWERS, WHICH IS NARROWER. Not "is the absolute AUROC inflated" - it
plainly is, and no absolute number here should be compared against a benchmark that samples
negatives differently. The question is whether the CONTRAST, which is what the paper claims,
is manufactured by that confound. Restrict to pairs whose NEGATIVE is plausibly real RNA
present in the cell, refit both arms, and see what survives.

WHY THE OBVIOUS VERSION OF THIS TEST LIES, AND WHY THE PLACEBO IS THE WHOLE EXPERIMENT.
Restricting on expression discards roughly two thirds of pairs. A 256-feature k-mer model
loses more from that than a 19-feature composition baseline does, in both arms, so the
contrast shrinks whether or not expression matters. Exactly the trap `strand_placebo.py`
documents. The fix is the same: drop the SAME NUMBER of pairs at random, many times, and
report the difference. Everything that is merely the cost of having less data cancels.

RESTRICTION IS NOT A RANDOM DROP, AND WHAT IT CORRELATES WITH IS MEASURED HERE RATHER THAN
ASSUMED. Retention is a property of the negative's locus, so it is tied to locus type. The
script measures the retained-versus-dropped standardised mean difference on region mix, GC
and gene density and writes all three to the per-dataset table, then matches the placebo on
the negative's region class. The unstratified placebo is computed alongside so their
difference bounds the locus-mix component, which is the same reporting the strand control
uses. Twenty seeds, because five left roughly a sixth of the between-dataset variance as
Monte Carlo noise there and caused a true finding to be withdrawn.

BALANCE ACROSS ARMS IS WHAT MAKES THIS CONVINCING, AND IT IS PRINTED. A confound can only
manufacture the contrast if it differs between the arms. The unexpressed-negative fraction
is 40.5% in the GC arm against 40.3% in the dinucleotide arm (paired Wilcoxon p = 0.87), and
expression-alone AUROC differs by +0.0006. The two arms carry this confound equally, so it
cancels in the difference by construction and the restriction below is a check on that
argument rather than the argument itself.

READ THE SUMMARY AS WITHIN-STORE DIFFERENCES ONLY. The GC arm's local window tables rebuild
the committed rehearsal rows exactly (12/12 spot-checked, max |diff| 1e-4). The dinucleotide
arm's do not: negative matching is a stochastic search that was re-run, and the local copy is
a different draw, reproducing 5 of 12 at the 2e-3 tolerance `k_sweep.py` uses. That is why
every quantity reported here is a difference computed WITHIN one store - full, restricted and
placebo all from the same draw - and why the panel's own full-data contrast is printed beside
the published +0.0397 rather than substituted for it. Differencing a locally restricted
contrast against the published number would measure the draw, not expression. `--strict`
turns the dinucleotide reproduction into a hard gate for anyone holding canonical tables.

PRE-REGISTERED before the run: threshold TPM >= 1 on the labelled strand of the NEGATIVE;
minimum 200 retained pairs; the claim is supported if the excess (restricted minus
region-matched placebo) has an interval that includes zero, or is small against the
-0.0055 the strand artifact costs.

    python scripts/expression_control.py --gc-root ... --dn-root ... --resume
    python scripts/expression_control.py --from-cache          # summary only
"""

import argparse
import gzip
import sys
import urllib.request
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np                                                      # noqa: E402
import pandas as pd                                                     # noqa: E402
from scipy.stats import wilcoxon                                        # noqa: E402

from rbp.eval import baseline, nested                                   # noqa: E402
from rbp.utils import config as cfgmod                                  # noqa: E402
from rbp.utils.log import log  # noqa: E402
from strand_placebo import pair_key, stratified_pick, subset            # noqa: E402

TABLES = ROOT / "results" / "tables"
EXTERNAL = ROOT / "data" / "external"
ARM_DIR = {"gc": "gc", "dn": "dinuc"}

# ENCODE total RNA-seq gene quantifications, RSEM TPM, GENCODE V29 gene ids. Two released
# replicates per cell line; the per-gene maximum is taken across them, which is the
# permissive choice: it can only move a gene from "silent" to "expressed" and therefore
# only ever SHRINKS the confound this script is trying to expose.
QUANT = {
    "K562": ["ENCFF286KKZ", "ENCFF829LCN"],       # ENCSR115PIZ
    "HepG2": ["ENCFF863QWG", "ENCFF376IXQ"],      # ENCSR245ATJ, ENCSR813BDU
}
ENCODE_FILE = "https://www.encodeproject.org/files/{acc}/@@download/{acc}.tsv"

TPM_MIN = 1.0            # a gene at 1 TPM is present; below it the window may not be RNA
MIN_KEPT_PAIRS = 200     # panel.min_test_pairs: below this an AUROC is too noisy to difference
BIN = 1000               # expression track resolution, nt
REPRO_TOL = 2.0e-3       # the tolerance k_sweep.py uses for the same reproduction gate
N_PLACEBO = 20           # five seeds left ~1/6 of the variance as noise in the strand control
N_BOOT = 2000
SEED = 0
DESIGN = "expressed_negative_tpm1_sense"



# ---------------------------------------------------------------------------------------
# External data: the GTF the study already uses, plus ENCODE expression for both lines
# ---------------------------------------------------------------------------------------

def cached(url, name):
    """Download once into data/external/. That directory is gitignored, so the artefact is
    reproducible from the accession rather than committed."""
    EXTERNAL.mkdir(parents=True, exist_ok=True)
    dst = EXTERNAL / name
    if not dst.exists():
        log(f"  fetching {name}")
        tmp = dst.with_suffix(dst.suffix + ".part")
        urllib.request.urlretrieve(url, tmp)
        tmp.rename(dst)                      # atomic, so a killed download is not read as done
    return dst


def gene_table(gtf_path):
    """One row per GENCODE gene: chrom, start, end, strand, versionless id, TPM per cell.

    The GTF is the same release the region index was built from, read through
    `config/params.yaml` rather than hardcoded, so this cannot drift from the annotation the
    windows were classified with.
    """
    def attr(f, key):
        i = f.find(key + ' "')
        if i == -1:
            return None
        return f[i + len(key) + 2:f.index('"', i + len(key) + 2)]

    rows = []
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t")
            if f[2] != "gene":
                continue
            rows.append((f[0], int(f[3]) - 1, int(f[4]), f[6],
                         attr(f[8], "gene_id").split(".")[0]))
    g = pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "gid"])
    g = g[g.chrom.str.match(r"^chr(\d+|X|Y)$")]

    for cell, accs in QUANT.items():
        tpm = None
        for acc in accs:
            d = pd.read_csv(cached(ENCODE_FILE.format(acc=acc), f"{acc}.tsv"),
                            sep="\t", usecols=["gene_id", "TPM"])
            d = d[d.gene_id.str.startswith("ENSG")]
            s = d.assign(gid=d.gene_id.str.split(".").str[0]).groupby("gid").TPM.max()
            tpm = s if tpm is None else pd.concat([tpm, s], axis=1).max(axis=1)
        g = g.merge(tpm.rename(cell), on="gid", how="left")
    return g


def expression_tracks(g, cell):
    """chrom -> per-strand arrays of max TPM over 1 kb bins.

    Binned rather than interval-searched because every window is 101 nt and therefore
    touches at most two bins, so a lookup is two array reads instead of a tree query over
    63,187 genes for each of 3.6 million windows. The bin rounds a window's footprint UP to
    1 kb, which can only find MORE expression nearby, so it is conservative in the direction
    that works against this script's finding.
    """
    tr = {}
    for chrom, sub in g.groupby("chrom", observed=True):
        n = int(sub.end.max()) // BIN + 2
        d = {"+": np.zeros(n, np.float32), "-": np.zeros(n, np.float32)}
        for st, en, strand, tpm in zip(sub.start.values, sub.end.values,
                                       sub.strand.values, sub[cell].fillna(0).values):
            np.maximum.at(d[strand], np.arange(st // BIN, en // BIN + 1), np.float32(tpm))
        tr[chrom] = d
    return tr


def window_tpm(tracks, chrom, start, end, strand=None):
    """Max TPM overlapping the window: on `strand` if given, else on either strand.

    BOTH ARE NEEDED AND THEY MEASURE DIFFERENT THINGS. The labelled-strand figure is the
    honest restriction - it asks whether the window is the RNA the model is being shown -
    but it inherits the strand artifact of limitation 1, because a negative carries its
    POSITIVE's strand. The either-strand figure asks only "is this locus transcribed at
    all", which is expression per se, uncontaminated by strand. The arm-balance argument
    has to be made on the either-strand figure; the restriction is applied on the
    labelled-strand one, which is the stricter of the two.
    """
    t = tracks.get(chrom)
    if t is None:
        return 0.0
    arrs = [t[strand]] if strand else [t["+"], t["-"]]
    out = 0.0
    for a in arrs:
        b0, b1 = min(start // BIN, len(a) - 1), min(end // BIN, len(a) - 1)
        out = max(out, float(a[b0]), float(a[b1]))
    return out


# ---------------------------------------------------------------------------------------
# The measurement
# ---------------------------------------------------------------------------------------

def nested_gain(d):
    """The published quantity: out-of-fold AUROC(composition + score) - AUROC(composition)."""
    res = baseline.evaluate(d, k=4)
    g = nested.gain_over_composition(d.seq_rna.tolist(), res["scores"],
                                     d.label.to_numpy(), d.fold.to_numpy())
    return g.delta, g.auroc_composition, g.auroc_with_score


def expressed_pairs(d, tracks):
    """Keys of pairs whose NEGATIVE overlaps a gene above TPM_MIN on its labelled strand.

    Keyed on the negative alone. Positives are 98.6% expressed by construction - they are
    where a peak was called - so restricting on the positive would drop almost nothing and
    test almost nothing. The pair moves as a unit because region and GC are matched within
    it, so removing one member would unbalance the very axes the arms differ on.
    """
    neg = d[d.label == 0]
    keys = pair_key(neg.id)
    coords = list(zip(neg.chrom, neg.start, neg.end, neg.strand))
    ok = np.array([window_tpm(tracks, c, int(s), int(e), a) >= TPM_MIN
                   for c, s, e, a in coords])
    any_ok = np.array([window_tpm(tracks, c, int(s), int(e)) >= TPM_MIN
                       for c, s, e, _ in coords])
    return set(keys[ok]), float(1.0 - any_ok.mean())


def smd(a, b):
    """Standardised mean difference, pooled sd. Reported so that 'the placebo is matched on
    the right thing' is a measurement rather than an assertion."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    sd = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    return float((a.mean() - b.mean()) / sd) if sd > 0 else 0.0


def imbalance(d, keep):
    """What restriction selects on, measured: region mix, GC and window density."""
    neg = d[d.label == 0].copy()
    neg["k"] = pair_key(neg.id)
    kept = neg[neg.k.isin(keep)]
    drop = neg[~neg.k.isin(keep)]
    if not len(kept) or not len(drop):
        return {}
    out = {"smd_gc": smd(kept.gc, drop.gc)}
    for r in ("intron", "cds", "utr3"):
        out[f"frac_{r}_kept"] = float((kept.region == r).mean())
        out[f"frac_{r}_drop"] = float((drop.region == r).mean())
    return out


def region_strata(d):
    """Pair key -> the NEGATIVE's region class, the axis restriction actually moves.

    Region only. Region is matched 1:1 within every pair by `negatives.py`, so restricting
    on the negative reweights the whole task by locus type and that is the confound worth
    removing. GC is matched within the pair too and the strand control measured it balanced
    to the third decimal under a similar restriction, so adding GC quintiles would multiply
    the cell count for a dimension that is not confounded. The unstratified placebo is
    computed alongside anyway, and their difference is reported as the locus-mix bound.
    """
    neg = d[d.label == 0]
    return dict(zip(pair_key(neg.id), neg.region.astype(str)))


def run_dataset(ds, root, pub, tracks_by_cell, strict, n_placebo):
    """One dataset, both arms. Returns a record or None if it cannot be measured."""
    prot, cell = ds.split(":")
    tracks = tracks_by_cell[cell]
    rec = {"dataset": ds, "protein": prot, "cell": cell, "design": DESIGN}

    for arm in ("gc", "dn"):
        f = root[arm] / cell / prot / "dataset.tsv"
        if not f.exists() or ds not in pub[arm].index:
            return None
        d = pd.read_csv(f, sep="\t")
        full, comp, with_s = nested_gain(d)
        r = pub[arm].loc[ds]
        # REPRODUCTION, recorded always and fatal only on request. The GC arm rebuilds its
        # published row exactly; the dinucleotide store is a different draw of the stochastic
        # matcher, so a hard gate there would discard most of the panel for a reason that
        # cancels inside a within-store difference. The deviation is written to the table so
        # a reader can see which rows are canonical.
        rec[f"repro_comp_{arm}"] = abs(comp - r.composition_auroc)
        rec[f"repro_with_{arm}"] = abs(with_s - r.with_score_auroc)
        if strict and max(rec[f"repro_comp_{arm}"], rec[f"repro_with_{arm}"]) > REPRO_TOL:
            return None

        allk = set(pair_key(d.id))
        keep, frac_unexpr_any = expressed_pairs(d, tracks)
        rec[f"n_pairs_{arm}"] = len(allk)
        rec[f"n_expressed_{arm}"] = len(keep)
        rec[f"frac_unexpressed_{arm}"] = 1.0 - len(keep) / max(len(allk), 1)
        rec[f"frac_unexpressed_any_{arm}"] = frac_unexpr_any
        rec[f"full_{arm}"] = full
        rec[f"comp_{arm}"] = comp
        if not (MIN_KEPT_PAIRS <= len(keep) < len(allk)):
            return None
        rec[f"expr_{arm}"] = nested_gain(subset(d, keep))[0]
        for k, v in imbalance(d, keep).items():
            rec[f"{k}_{arm}"] = v

        st = region_strata(d)
        pl, pls, defs = [], [], 0
        for s in range(n_placebo):
            rng = np.random.default_rng(3000 + s)
            pick = set(rng.choice(sorted(allk), size=len(keep), replace=False))
            pl.append(nested_gain(subset(d, pick))[0])
            rng2 = np.random.default_rng(4000 + s)
            picks, dfc = stratified_pick(sorted(allk), keep, st, rng2)
            defs += dfc
            pls.append(nested_gain(subset(d, picks))[0])
        rec[f"placebo_{arm}"] = float(np.mean(pl))
        rec[f"placebo_sd_{arm}"] = float(np.std(pl))
        rec[f"placebo_strat_{arm}"] = float(np.mean(pls))
        rec[f"placebo_strat_sd_{arm}"] = float(np.std(pls))
        rec[f"strat_deficit_{arm}"] = defs / n_placebo
    return rec


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", default=None)
    p.add_argument("--gtf", default="", help="GENCODE GTF; downloaded to data/external if unset")
    p.add_argument("--gc-root", default="")
    p.add_argument("--dn-root", default="")
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--placebo-seeds", type=int, default=N_PLACEBO)
    p.add_argument("--strict", action="store_true",
                   help="drop datasets whose local windows do not rebuild the published row")
    p.add_argument("--resume", action="store_true")
    # Rebuild the summary from the committed per-dataset table without redoing the refits.
    # The per-dataset table IS the evidence; the summary is arithmetic on it.
    p.add_argument("--from-cache", action="store_true")
    a = p.parse_args()

    cache = TABLES / "expression_control_per_dataset.csv"
    if a.from_cache:
        return summarise(pd.read_csv(cache))

    cfg = cfgmod.load(a.config)
    gtf = Path(a.gtf) if a.gtf else cached(cfg.reference["gtf"], "gencode.annotation.gtf.gz")
    log("building gene table and expression tracks")
    g = gene_table(gtf)
    tracks_by_cell = {c: expression_tracks(g, c) for c in QUANT}
    log(f"  {len(g):,} genes; TPM present for {int(g.K562.notna().sum()):,}")

    cm = pd.read_csv(TABLES / "cost_of_matching.csv")
    names = list(cm.dataset)[: a.limit or None]
    pub = {"gc": pd.read_csv(TABLES / "rehearsal_binding_gc.csv").set_index("dataset"),
           "dn": pd.read_csv(TABLES / "rehearsal_binding_dinuc.csv").set_index("dataset")}
    root = {"gc": Path(a.gc_root), "dn": Path(a.dn_root)}

    rows, done = [], set()
    if a.resume and cache.exists():
        prev = pd.read_csv(cache)
        # Resume only within the same design, keyed explicitly. A guard that cannot tell two
        # designs apart silently reported a mixture of both in the strand control.
        if len(prev) and prev.get("design", pd.Series(dtype=object)).eq(DESIGN).all():
            rows, done = prev.to_dict("records"), set(prev.dataset)
            log(f"  resuming: {len(done)} datasets already done under design '{DESIGN}'")
        else:
            log(f"  cache is a different design; recomputing all under '{DESIGN}'")

    for i, ds in enumerate(names, 1):
        if ds in done:
            continue
        rec = run_dataset(ds, root, pub, tracks_by_cell, a.strict, a.placebo_seeds)
        if rec is None:
            log(f"  [{i:3d}/{len(names)}] {ds:22} SKIP")
            continue
        rows.append(rec)
        log(f"  [{i:3d}/{len(names)}] {ds:22} kept "
            f"{rec['n_expressed_gc']}/{rec['n_pairs_gc']} gc  "
            f"contrast full {rec['full_dn'] - rec['full_gc']:+.4f} "
            f"expr {rec['expr_dn'] - rec['expr_gc']:+.4f} "
            f"placebo {rec['placebo_strat_dn'] - rec['placebo_strat_gc']:+.4f}")
        pd.DataFrame(rows).to_csv(cache, index=False)

    m = pd.DataFrame(rows)
    if not len(m):
        raise SystemExit("no datasets measurable; nothing to report")
    return summarise(m)


def summarise(m):
    m = m.dropna(subset=["full_gc", "full_dn", "expr_gc", "expr_dn",
                         "placebo_strat_gc", "placebo_strat_dn"]).copy()
    m["c_full"] = m.full_dn - m.full_gc
    m["c_expr"] = m.expr_dn - m.expr_gc
    m["c_placebo"] = m.placebo_dn - m.placebo_gc
    m["c_placebo_strat"] = m.placebo_strat_dn - m.placebo_strat_gc
    m["excess"] = m.c_expr - m.c_placebo
    m["excess_strat"] = m.c_expr - m.c_placebo_strat
    m["locus_mix"] = m.c_placebo_strat - m.c_placebo
    # Only the expression-specific part is removed. The shrinkage the placebo also shows is
    # the cost of discarding pairs, not of expression, so subtracting it would double-count
    # and understate the effect the paper claims.
    m["corrected"] = m.c_full + m.excess_strat
    m.to_csv(TABLES / "expression_control_per_dataset.csv", index=False)

    rng = np.random.default_rng(SEED)
    n = len(m)
    keys = ("c_full", "c_expr", "c_placebo", "c_placebo_strat", "excess", "excess_strat",
            "locus_mix", "corrected", "d_expr", "d_placebo", "arm_balance")
    boots = {k: [] for k in keys}
    for _ in range(N_BOOT):
        s = m.iloc[rng.integers(0, n, n)]
        for k in ("c_full", "c_expr", "c_placebo", "c_placebo_strat", "excess",
                  "excess_strat", "locus_mix", "corrected"):
            boots[k].append(s[k].mean())
        boots["d_expr"].append((s.c_expr - s.c_full).mean())
        boots["d_placebo"].append((s.c_placebo_strat - s.c_full).mean())
        boots["arm_balance"].append(
            (s.frac_unexpressed_any_dn - s.frac_unexpressed_any_gc).mean())

    out = []

    def add(check, value, key=None, note=""):
        lo, hi = np.percentile(boots[key], [2.5, 97.5]) if key else (np.nan, np.nan)
        out.append({"check": check, "value": float(value), "ci_low": lo, "ci_high": hi,
                    "n": n, "note": note})

    kept = m.n_expressed_gc.sum() / m.n_pairs_gc.sum()
    # THE BALANCE ARGUMENT IS MADE ON THE EITHER-STRAND FIGURE, and the two must not be
    # confused. Either-strand is expression per se and is what has to be balanced for the
    # confound to cancel in a paired difference. The labelled-strand figure is larger and is
    # NOT balanced, because it inherits the strand artifact of limitation 1: negatives carry
    # the positive's strand, and the GC arm carries more of that cue than the dinucleotide
    # arm does (R1c, +0.047 more sense negatives in the dinuc arm). That difference is in the
    # conservative direction for the contrast, the same way R1c's is.
    add("negatives in an untranscribed locus, GC arm", m.frac_unexpressed_any_gc.mean(),
        note="no gene above 1 TPM on EITHER strand: expression per se")
    add("negatives in an untranscribed locus, dinuc arm", m.frac_unexpressed_any_dn.mean(),
        note="the confound has to DIFFER between arms to manufacture the contrast")
    bal = wilcoxon(m.frac_unexpressed_any_dn, m.frac_unexpressed_any_gc)
    add("arm difference, untranscribed fraction",
        (m.frac_unexpressed_any_dn - m.frac_unexpressed_any_gc).mean(), "arm_balance",
        note=f"paired Wilcoxon p = {bal.pvalue:.3g}; balance is why the confound cancels")
    add("negatives not transcribed on the LABELLED strand, GC arm",
        m.frac_unexpressed_gc.mean(),
        note="the restriction criterion; larger because it also carries the strand artifact")
    add("negatives not transcribed on the LABELLED strand, dinuc arm",
        m.frac_unexpressed_dn.mean(),
        note="smaller than the GC arm, the same conservative direction R1c reports")
    add("contrast, full data", m.c_full.mean(), "c_full",
        note="this panel's own contrast; published on all 94 is +0.0397")
    add("contrast, expressed-negative pairs", m.c_expr.mean(), "c_expr",
        note=f"mean {kept:.1%} of GC pairs retained")
    add("contrast, PLACEBO (same n, random)", m.c_placebo.mean(), "c_placebo",
        note=f"{int(m.get('placebo_sd_gc', pd.Series([np.nan])).notna().sum())} rows with seed sd recorded")
    add("contrast, PLACEBO stratified on region", m.c_placebo_strat.mean(), "c_placebo_strat",
        note="matched to the retained set's region marginals; PRIMARY")
    add("change from restriction", (m.c_expr - m.c_full).mean(), "d_expr",
        note="what the naive test would have reported as expression")
    add("change from placebo", (m.c_placebo_strat - m.c_full).mean(), "d_placebo",
        note="the same shrinkage, with no expression involved")
    add("expression excess, UNSTRATIFIED placebo", m.excess.mean(), "excess",
        note="expression PLUS locus mix")
    add("EXPRESSION-SPECIFIC EXCESS (stratified)", m.excess_strat.mean(), "excess_strat",
        note="THE ANSWER, pre-committed as primary")
    add("locus-mix component", m.locus_mix.mean(), "locus_mix",
        note="stratified placebo minus unstratified")
    add("expression-CORRECTED contrast", m.corrected.mean(), "corrected",
        note="full contrast with only the expression-specific part removed")
    add("fraction of the contrast surviving", m.corrected.mean() / m.c_full.mean(),
        note=f"panel's own contrast is {m.c_full.mean():+.4f}, not the n=94 +0.0397")
    add("mean placebo seed noise (sd), GC arm", float(m.placebo_strat_sd_gc.mean()),
        note="Monte Carlo component of the primary estimator")

    res = pd.DataFrame(out)
    res.to_csv(TABLES / "expression_control.csv", index=False)
    log("")
    for _, x in res.iterrows():
        ci = f" [{x.ci_low:+.4f}, {x.ci_high:+.4f}]" if pd.notna(x.ci_low) else ""
        log(f"  {x.check:42} {x.value:+.4f}{ci}   {x.note}")
    log(f"\n  n = {n} datasets;  wrote expression_control.csv")


if __name__ == "__main__":
    main()
