"""A source line ending mid-word in a hyphen becomes two words, and the PDF may not show it.

WHY. `protocol-` on one line and `specific` on the next is not `protocol-specific`: LaTeX turns
the newline into an interword space and sets `protocol- specific`. It has happened twice here,
in results.tex and again in introduction.tex, and neither was caught by reading the PDF. The
reason is that TeX may break the line exactly at that space, in which case the page shows
`protocol-` at the right margin and `specific` at the left, which is indistinguishable from
correct hyphenation. The defect is then invisible until the paragraph reflows, and it is always
wrong in the text layer, so copy-paste and full-text search see the space.

So this checks the SOURCE, where the property is decidable, rather than the render, where it is
not. A hyphen ending a line whose next line continues in lower case is the signature; a hyphen
after a space is a minus sign or a dash and is left alone.
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAN = ROOT / "manuscript"

# A word character immediately before the hyphen. `... score) -` at the end of a display
# equation in Methods is a minus sign, not a broken word, and has a space in front of it.
BROKEN = re.compile(r"[A-Za-z0-9]-$")


def offenders(path):
    lines = path.read_text().split("\n")
    out = []
    for i, line in enumerate(lines):
        if line.lstrip().startswith("%") or not BROKEN.search(line.rstrip()):
            continue
        nxt = next((x for x in lines[i + 1:] if x.strip()), "")
        # A following line that starts a new sentence, a macro or a comment is a real dash at a
        # line end; only a lower-case continuation is a word that was split.
        if re.match(r"[a-z]", nxt.lstrip()):
            out.append((i + 1, line.rstrip(), nxt.strip()[:40]))
    return out


def main():
    srcs = sorted(MAN.glob("*.tex")) + sorted((MAN / "sections").glob("*.tex"))
    if not srcs:
        sys.exit("no manuscript sources found; this check looked at NOTHING")
    bad = [(p, o) for p in srcs for o in [offenders(p)] if o]
    for p, os_ in bad:
        for ln, line, nxt in os_:
            print(f"{p.relative_to(ROOT)}:{ln}: word split across lines by a hyphen\n"
                  f"    ...{line[-58:]}\n    {nxt}...")
    n = sum(len(o) for _, o in bad)
    if n:
        print(f"\n{n} split word(s). Join the hyphen to the word that follows it.")
        return 1
    print(f"  {len(srcs)} manuscript sources, no word split across a line by a hyphen")
    return 0


if __name__ == "__main__":
    sys.exit(main())
