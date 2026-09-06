"""Which datasets are in the study, for a given negative arm.

WHY THIS EXISTS AS A SHARED MODULE. Three analysis scripts each had their own copy of
"read the panel file, filter on pairs, build a path to dataset.tsv". All three read
`panel_final_<cell>.tsv` and all three then loaded from `data/processed/` -- the GC
directory. That was wrong in two directions at once once the dinucleotide arm existed:

  * the filename is not keyed by arm, so whichever arm ran last overwrote the other's
    pair counts. Measured: AATF's file said 1081 (its dinuc count) while
    its GC dataset holds 1043 pairs.
  * the dataset directory was hardcoded to the GC arm, so the dinuc counts were used to
    filter GC data.

prepare.py already carried the right reasoning for the other axis -- "keyed by cell line:
a single panel_final.tsv would have the second cell line overwrite the first and silently
halve the study" -- and it simply was not carried across to the arm. One implementation,
used everywhere, is what stops that happening a third time.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]

# The arm decides BOTH the panel file and the directory the datasets live in. Keeping the
# pair together in one place is the point: they went out of step precisely because two
# different lines of code chose them independently.
ARMS = {
    "gc": "data/processed",
    "dinuc": "data/processed_dinucmatch",
    # The bias-aware arm. Its windows are built by scripts/build_neg2.py directly under
    # processed/neg2; this entry is what lets --arm neg2 resolve through the path helpers.
    "neg2": "data/processed_neg2",
    # The REGION-MATCHED bias-aware arm, built by build_neg2.py --match-region. Same 94
    # datasets and the same 456,734 pairs as `neg2`, drawn stratified on transcript region
    # instead of uniformly inside the fold. It exists because `neg2` matches fold only, and
    # region alone separates that arm's classes at median AUROC 0.748. It is a separate arm
    # rather than a replacement: the published result is the unstratified draw, and this one
    # bounds how much of it is region mix.
    "neg2_rm": "data/processed_neg2_rm",
}

# The two composition-matched arms, which is what "both arms" means everywhere in this project
# that predates neg2. Stages written against that pairing should say so explicitly rather than
# iterate ARMS and silently acquire a third arm they were never designed for.
COMPOSITION_MATCHED_ARMS = ("dinuc", "gc")


def check_arm(arm):
    if arm not in ARMS:
        raise ValueError(f"unknown negative arm {arm!r}; expected one of {sorted(ARMS)}")
    return arm


def panel_path(cell, arm, root=None):
    """config/panel_final_<cell>_<arm>.tsv"""
    check_arm(arm)
    return (Path(root) if root else ROOT) / "config" / f"panel_final_{cell}_{arm}.tsv"


def excluded_path(cell, arm, root=None):
    check_arm(arm)
    return (Path(root) if root else ROOT) / "config" / f"panel_excluded_{cell}_{arm}.tsv"


def data_dir(cell, arm, root=None):
    check_arm(arm)
    return (Path(root) if root else ROOT) / ARMS[arm] / cell


def read_panel(cell, arm, root=None):
    """[(protein, pairs)] from the frozen panel file, or [] if it has not been written.

    Parsed by hand rather than with pandas so that `rbp.utils` stays importable in the
    CPU container without pulling pandas into a path that does not otherwise need it.
    """
    f = panel_path(cell, arm, root)
    if not f.exists():
        return []
    out = []
    for ln in f.read_text().strip().splitlines()[1:]:
        parts = ln.split("\t")
        out.append((parts[0], int(parts[2])))
    return out


def datasets(cells, arm, min_pairs=None, root=None, require_files=True):
    """'PROTEIN:CELL' -> path to dataset.tsv, for everything in the arm's final panel.

    A protein listed in the panel whose dataset file is missing is skipped rather than
    raising, because the panel is a record of what qualified and the files can legitimately
    have been cleaned up. The count is returned alongside so a caller can notice.

    `require_files=False` returns the intended set regardless of what is on this disk. The
    container has config/ but not data/, and the cloud driver builds its manifest from the
    panel alone; without this the same call would answer "the study is empty" there and
    "the study is 189 datasets" on the laptop.
    """
    out, missing = {}, []
    for cell in cells:
        for protein, pairs in read_panel(cell, arm, root):
            if min_pairs is not None and pairs < min_pairs:
                continue
            p = data_dir(cell, arm, root) / protein / "dataset.tsv"
            if p.exists() or not require_files:
                out.setdefault(f"{protein}:{cell}", p)
            if not p.exists():
                missing.append(f"{protein}:{cell}")
    return out, missing


def arm_of(cfg, arm=None):
    """The negative arm to use: an explicit override, else the config's primary."""
    return check_arm(arm or cfg["negatives"]["primary_arm"])


def cells_of(cfg):
    return list(cfg.encode["cell_lines"])


def study(cfg, arm=None, root=None, require_files=True):
    """THE definition of what the study runs on. Everything else asks this.

    Returns ({'PROTEIN:CELL': path_to_dataset_tsv}, [missing keys]).

    Four scripts used to answer this question independently, and all four got it wrong in
    the same way: they read `config/panel_final.tsv`, a file keyed by neither cell line nor
    arm, left over from a 17-protein development panel. A sweep driven by it would have
    trained on 17 proteins in one cell line on the GC arm, and reported it as the study.

    Cells, arm and the pair threshold all come from config, so changing the study is a
    config edit and cannot be done by half.
    """
    return datasets(cells_of(cfg), arm_of(cfg, arm), cfg.cv["min_pairs"], root,
                    require_files)


def write_panel(cell, arm, keep, drop, root=None):
    """Write both panel files for one (cell, arm). `keep` is [(protein, pairs)];
    `drop` is [(protein, pairs, reason)]."""
    check_arm(arm)
    pf, xf = panel_path(cell, arm, root), excluded_path(cell, arm, root)
    pf.parent.mkdir(parents=True, exist_ok=True)
    pf.write_text("protein\tcell_line\tpairs\n"
                  + "".join(f"{p}\t{cell}\t{n}\n" for p, n in sorted(keep)))
    xf.write_text("protein\tcell_line\tpairs\treason\n"
                  + "".join(f"{p}\t{cell}\t{n}\t{w}\n" for p, n, w in sorted(drop)))
    return pf, xf
