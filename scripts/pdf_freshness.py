"""The tracked PDFs must be the output of the tracked sources.

    python scripts/pdf_freshness.py

WHY THIS EXISTS. CI proves that LaTeX can build the sources. It has never proved that the
binary a reader downloads is what those sources produce. Those are different claims, and on
2026-09-07 they came apart: an external audit built the audited commit from scratch and found
the fresh text differing from `manuscript/paper.pdf` on nine pages.

The audit read that as a stale PDF and it was worse than that. Two commits had bumped a stated
page count with a blind substitution, `57` to `58` and then `58` to `59`, and in a regex
`\\b58\\b` matches the `58` inside `6.58` because `.` is not a word character. Fourteen unrelated
scientific numbers moved. The committed PDF still held the CORRECT ones. So the tracked binary
and the tracked source disagreed, and the source was the wrong half.

Either direction of that disagreement is a research-object integrity failure, and nothing in
this repository could see it. This script is the check that can.

How it compares, and why not bytes. pdflatex embeds a build timestamp and a document ID, so two
builds of identical sources differ in bytes. Demanding byte identity would produce a gate that
fails for a reason nobody should act on. Instead the PDF is compared on what a reader actually
reads:

  * the full sequence of NUMERIC TOKENS, in order, which is the thing that went wrong and the
    thing that matters. A single moved digit fails.
  * the whitespace-normalised text, reported as a similarity ratio, with the first differing
    span printed. This catches reworded prose that carries no number.
  * the page count.

The build happens in a temporary copy, so the tree is never modified and a failure leaves
nothing to clean up.
"""

import difflib
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.utils.log import log  # noqa: E402

MANUSCRIPT = ROOT / "manuscript"
DOCS = ("paper.pdf", "supplementary.pdf")


def text_of(pdf):
    from pypdf import PdfReader
    r = PdfReader(str(pdf))
    return len(r.pages), "".join(p.extract_text() or "" for p in r.pages)


def norm(t):
    return re.sub(r"\s+", " ", t).strip()


def numbers(t):
    return re.findall(r"\d+(?:\.\d+)?", t)


def main():
    for d in DOCS:
        if not (MANUSCRIPT / d).exists():
            sys.exit(f"{d} is not committed, so there is nothing to check it against")

    tmp = Path(tempfile.mkdtemp(prefix="pdf-fresh-"))
    bad = []
    try:
        # build.sh reads figures as ../results/figures/, so the temp tree has to mirror the
        # repository layout, not just hold a copy of manuscript/. Figures are symlinked rather
        # than copied: they are large, and this script must not be able to modify them.
        build = tmp / "manuscript"
        shutil.copytree(MANUSCRIPT, build)
        (tmp / "results").mkdir(exist_ok=True)
        for name in ("figures", "tables"):
            src = ROOT / "results" / name
            if src.exists():
                (tmp / "results" / name).symlink_to(src)
        for d in DOCS:                      # so a failed build cannot pass by leaving the old file
            (build / d).unlink()
        r = subprocess.run(["./build.sh"], cwd=str(build), capture_output=True, text=True)
        if r.returncode != 0:
            tail = (r.stderr or r.stdout).strip().splitlines()[-6:]
            sys.exit("the clean build FAILED:\n  " + "\n  ".join(tail))

        for d in DOCS:
            if not (build / d).exists():
                bad.append(f"{d}: the clean build did not produce it")
                continue
            n_have, t_have = text_of(MANUSCRIPT / d)
            n_fresh, t_fresh = text_of(build / d)

            if n_have != n_fresh:
                bad.append(f"{d}: tracked has {n_have} pages, a clean build gives {n_fresh}")

            a, b = numbers(t_have), numbers(t_fresh)
            if a != b:
                sm = difflib.SequenceMatcher(None, a, b)
                shown = 0
                for tag, i1, i2, j1, j2 in sm.get_opcodes():
                    if tag == "equal" or shown >= 6:
                        continue
                    bad.append(f"{d}: NUMBERS {tag} at token {i1}: "
                               f"{a[i1:i2][:4]} -> {b[j1:j2][:4]}")
                    shown += 1

            ha, hb = norm(t_have), norm(t_fresh)
            if ha != hb:
                ratio = difflib.SequenceMatcher(None, ha, hb).quick_ratio()
                i = next((k for k in range(min(len(ha), len(hb))) if ha[k] != hb[k]),
                         min(len(ha), len(hb)))
                bad.append(f"{d}: TEXT differs (similarity {ratio:.4f}); first at char {i}: "
                           f"{ha[i:i + 60]!r} -> {hb[i:i + 60]!r}")
            else:
                log(f"  {d}: {n_have} pages, {len(a)} numeric tokens, identical to a clean build")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    if bad:
        log("\n  THE TRACKED PDF IS NOT THE OUTPUT OF THE TRACKED SOURCE:\n")
        for m in bad[:14]:
            log(f"    {m}")
        if len(bad) > 14:
            log(f"    ... and {len(bad) - 14} more")
        log("\n  Rebuild with `cd manuscript && ./build.sh` and commit the result. If the "
            "SOURCE is what changed by accident, fix the source first: a fresh build is not "
            "automatically the correct half.")
        sys.exit(1)
    log("  both tracked PDFs are the output of the tracked sources")


if __name__ == "__main__":
    main()
