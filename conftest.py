"""Make the tests import THIS repository's `rbp`, not the sibling project's.

The venv is shared with `../rna-binding-proteins`, and its editable install writes that
project's `src` onto sys.path via a .pth file. Most test modules open with

    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

but several -- test_annotation.py among them -- do not, and pytest imports modules in
alphabetical order. So `rbp` was bound to whichever tree the FIRST module to import it
happened to reach, which was the sibling project's, and every later insert was too late
because the package was already in sys.modules.

The two trees are near-identical, which is why this went unnoticed: they differ in
`utils/cloud.py`, `variants/conservation.py` and `variants/phylop.py`. That is not a reason
to leave it. A suite that silently tests a different checkout than the one being committed
cannot support any claim made about this repository, and the divergence is only ever going
to grow.

A rootdir conftest is loaded before any test module, so putting this repo's src first here
settles it once for every test regardless of collection order.
"""

import sys
from pathlib import Path

SRC = str(Path(__file__).resolve().parent / "src")
if sys.path[:1] != [SRC]:
    sys.path.insert(0, SRC)
