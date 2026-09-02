"""Rebuild the sweep's input store on local disk, so a run needs no cloud bucket at all.

    python scripts/build_store.py --store ../rbp-store --arm gc

WHY. `gs://rbp-repro-2026-derived` is unreadable: the billing account was closed. Everything
the GPU sweep reads is still here -- the frozen panels are in `config/`, and the prepared
windows are in the sibling project that produced them -- so the store can be reconstituted
rather than recovered.

Datasets are SYMLINKED, not copied. The GC arm is 244 MB and this laptop has under 10 GB
free; the store is a naming layer over files that already exist, not a second copy of them.

Writes `panel/<arm>/panel_final_<cell>_<arm>.tsv` and
`processed/<arm>/<cell>/<protein>/dataset.tsv`, which is exactly what cloud_train.py's
`manifest` and `run` modes read.
"""

import argparse
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.utils import panel as panelmod  # noqa: E402

CELLS = ("K562", "HepG2")
# Where the prepared windows actually live. The rebuild ran in this repo but wrote its
# processed arms into the project it grew out of, and those directories are the originals:
# AATF:K562 gc is 2086 rows, matching config/panel_final_K562_gc.tsv (1043 pairs) and
# rehearsal_binding_gc.csv (n=2086).
DEFAULT_DATA_ROOT = ROOT.parent / "rna-binding-proteins"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", required=True)
    p.add_argument("--arm", default="gc", choices=sorted(panelmod.ARMS))
    p.add_argument("--data-root", default=str(DEFAULT_DATA_ROOT))
    a = p.parse_args()

    store, data_root = Path(a.store), Path(a.data_root)
    arm_dir = data_root / panelmod.ARMS[a.arm]
    if not arm_dir.exists():
        sys.exit(f"no prepared data at {arm_dir}")

    n_ds = 0
    missing = []
    for cell in CELLS:
        src = ROOT / "config" / f"panel_final_{cell}_{a.arm}.tsv"
        if not src.exists():
            sys.exit(f"no frozen panel at {src}")
        dst = store / "panel" / a.arm / src.name
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)

        for protein, _pairs in panelmod.read_panel(cell, a.arm, root=ROOT):
            ds = arm_dir / cell / protein / "dataset.tsv"
            if not ds.exists():
                missing.append(f"{protein}:{cell}")
                continue
            link = store / "processed" / a.arm / cell / protein / "dataset.tsv"
            link.parent.mkdir(parents=True, exist_ok=True)
            if link.is_symlink() or link.exists():
                link.unlink()
            link.symlink_to(ds.resolve())
            n_ds += 1

    print(f"store {store}  arm={a.arm}")
    print(f"  {n_ds} datasets linked from {arm_dir}")
    if missing:
        print(f"  {len(missing)} in the panel but not prepared: {missing[:8]}"
              f"{' ...' if len(missing) > 8 else ''}")


if __name__ == "__main__":
    main()
