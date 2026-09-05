"""Counts a release document states about the release must be derived from the release.

    python scripts/release_consistency.py

WHY THIS EXISTS. An external review read the repository and found eight stale figures in the
submission package: 28 pages against a 48-page PDF, 20 references against 26, six main figures
against seven, Tables 1 to 7 against fourteen, 768 numeric assertions against 937, 696 in the
Zenodo template, and a recomputation error of 2.2e-16 against the manuscript's 3.3e-16. Every
one had been correct when it was typed. Every one was copied by hand into a file that no gate
read.

WHY THE EXISTING GATES DID NOT CATCH THEM. scripts/audit_manuscript.py scans prose for numbers
with no source, and it now scans these documents too, but it cannot catch these: 28, 20 and 7
are small integers that collide with something in a 855-value haystack, and 2.2e-16 is written
in an exponent form its regex does not match. Orphan-hunting asks "is this number real
somewhere". The question here is the different and stricter one: "is this number still true of
THIS artefact". That needs the artefact, not a haystack.

SO THE FACTS ARE DERIVED AND THE PROSE IS MATCHED AGAINST THEM. Each fact below is computed
from the built PDF, the LaTeX source, or a committed table. Each pattern says where a document
is allowed to state that fact. A mismatch is an error and a fact that no document states is
reported too, because a fact nobody quotes is a pattern that has silently stopped matching --
which is how a check like this dies quietly rather than failing.
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.utils.log import log  # noqa: E402

MANUSCRIPT = ROOT / "manuscript"
PDF = MANUSCRIPT / "paper.pdf"
SECTIONS = sorted((MANUSCRIPT / "sections").glob("*.tex"))
TEX = [MANUSCRIPT / "paper.tex"] + SECTIONS

# The documents that make claims about the release. The run chronicle under docs/ is excluded
# for the reason given in audit_manuscript.py: it is pasted terminal output, and every integer
# in it is an HTTP status or a task index.
DOCS = [ROOT / "README.md", ROOT / "SUBMISSION.md",
        ROOT / "docs" / "REPRODUCE.md", ROOT / "docs" / "PANELS.md",
        ROOT / "docs" / "ZENODO.md"] + TEX


def _tex(paths=None):
    return "\n".join(p.read_text() for p in (paths or TEX))


def pdf_pages():
    try:
        import pypdf
    except ImportError:
        return None
    return len(pypdf.PdfReader(str(PDF)).pages) if PDF.exists() else None


def n_references():
    b = MANUSCRIPT / "sections" / "bibliography.tex"
    if not b.exists():
        return None
    # \bibitem is the manual form; a thebibliography built by hand uses it exclusively here.
    return len(re.findall(r"\\bibitem", b.read_text()))


def n_environments(kind):
    """figure/table environments, starred or not, in the typeset source."""
    return len(re.findall(rf"\\begin\{{{kind}\*?\}}", _tex()))


def n_figure_files():
    d = MANUSCRIPT / "figures"
    return len(list(d.glob("*.pdf"))) if d.exists() else None


def n_verify_checks():
    """The verifier's own total, read from the summary it writes.

    Read from the artefact and not from verify.py's source, because the count is a property
    of a run: gates skip when a table is absent, and a number taken from the source would be
    the number the code COULD reach rather than the one it did.
    """
    f = ROOT / "results" / "tables" / "verify_summary.csv"
    if not f.exists():
        return None
    for ln in f.read_text().strip().splitlines()[1:]:
        k, _, v = ln.partition(",")
        if k.strip() == "assertions":
            return int(v)
    return None


def abstract_words():
    """Words in the abstract, counting a maths token as the one word a form's counter sees.

    Counted the way the bioRxiv textarea counts and not the way `wc -w` does: control
    sequences vanish, `$+0.0137$` is one word rather than three, and a brace group contributes
    only its text. A naive strip returned 618 for a 350-word abstract, which would have failed
    the release on the checker's arithmetic rather than on the document's.
    """
    m = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}",
                  (MANUSCRIPT / "paper.tex").read_text(), re.S)
    if not m:
        return None
    body = re.sub(r"\$[^$]*\$", " X ", m.group(1))      # one token per maths span
    body = re.sub(r"\\[a-zA-Z]+\*?", " ", body)         # control sequences are not words
    body = re.sub(r"[{}~\\]", " ", body)
    return len([w for w in body.split() if re.search(r"[A-Za-z0-9]", w)])


def collected_tests():
    """Tests pytest actually collects here, not a number typed into the README.

    THIS COUNT MOVES WHEN A SCRIPT IS ADDED, not only when a test is. Two suites parametrise
    over files -- test_no_hardcoded_project.py over every tracked .py, .sh, .yaml and .json,
    test_figure_output.py over every figure PDF -- so adding scripts/common_positives.py adds a
    test case without anyone writing one. That is surprising the first time and is why the
    number is derived here rather than trusted where it is typed.

    Returns None where the whole suite cannot be collected -- no torch, or a caller passing
    --ignore -- for the reason tests/unit/test_suite_size.py gives at length: a subset compared
    against a whole-suite figure is a meaningless comparison that fails builds for the wrong
    reason. golden.yaml's `min_tests_passing` is the floor; this is the census.
    """
    import subprocess
    r = subprocess.run([sys.executable, "-m", "pytest", "tests", "--collect-only", "-q"],
                       cwd=ROOT, capture_output=True, text=True,
                       env={**__import__("os").environ, "PYTHONPATH": str(ROOT / "src"),
                            "GOOGLE_CLOUD_PROJECT": "release-consistency-no-such-project"})
    # SUMMED FROM THE PER-FILE LINES, not read off a summary line. pyproject sets addopts="-q"
    # already, so the -q here makes it -qq and pytest drops the "N tests collected" line
    # entirely, printing only `path: count`. Matching the summary silently returned None on
    # every run -- a skip that looks like a missing artefact rather than a broken parse, which
    # is the failure mode this whole script exists to catch.
    #
    # Gated on pytest's exit status and not on the word "error" appearing in stdout: this
    # suite has test ids containing "error".
    counts = [int(m) for m in re.findall(r"^tests/.*?: (\d+)$", r.stdout, re.M)]
    return sum(counts) if counts and r.returncode == 0 else None


FACTS = {
    "tests collected": (collected_tests, [
        r"#\s*(\d+),? (?:passed|needs torch)",
        r"tests/\s+(\d+) tests",
    ]),
    "pages": (pdf_pages, [
        r"paper\.pdf`? \((\d+) pages\)",
        r"(\d+)[- ]page (?:manuscript|PDF|paper)",
    ]),
    "references": (n_references, [
        r"then the declarations and (\d+) references",
        r"(\d+) references\b",
    ]),
    "figure environments": (lambda: n_environments("figure"), [
        r"(\d+) (?:main )?figures? are typeset",
        r"the (\d+) main figures",
    ]),
    "table environments": (lambda: n_environments("table"), [
        r"Tables 1 to (\d+) are typeset",
    ]),
    "verify checks": (n_verify_checks, [
        r"(\d{3,4})\s+(?:numeric\s+)?(?:verification\s+)?assertions",
        r"All (\d{3,4}) verification checks",
        r"(\d{3,4})/\d{3,4} checks",
    ]),
    "abstract words": (abstract_words, [
        r"(\d+) words, no markup",
    ]),
}

# A word count is not a checksum: LaTeX macro stripping and a form's own counter will differ by
# a word or two, and failing a release on that would train everyone to ignore this script.
TOLERANCE = {"abstract words": 4}


def main():
    problems, unstated, facts = [], [], []
    for name, (derive, patterns) in FACTS.items():
        truth = derive()
        if truth is None:
            log(f"  {name:22} SKIP (artefact absent)")
            continue
        facts.append((name, truth))
        tol = TOLERANCE.get(name, 0)
        seen = 0
        for doc in DOCS:
            if not doc.exists():
                continue
            for i, line in enumerate(doc.read_text().splitlines(), 1):
                # Deduplicated on where the number sits, not on which pattern found it. Two
                # patterns for the same fact can both match one line -- a general one and the
                # specific phrasing it generalises -- and reporting that as two stale claims
                # doubles the count of a single edit.
                for start in sorted({m.start(1) for pat in patterns
                                     for m in re.finditer(pat, line)}):
                    seen += 1
                    got = int(re.match(r"\d+", line[start:]).group(0))
                    if abs(got - truth) > tol:
                        problems.append(
                            f"{doc.relative_to(ROOT)}:{i} says {name} = {got}, "
                            f"derived from the release: {truth}\n      {line.strip()[:110]}")
        log(f"  {name:22} {truth:<8} stated in {seen} place(s)")
        if not seen:
            unstated.append(name)

    # THE FACTS ARE COMMITTED, which is what makes them quotable. scripts/audit_manuscript.py
    # flags any number in a released document that no table can source, and a page count or a
    # test census lives in no result table -- so widening that audit to README and SUBMISSION
    # correctly reported them as unsourced. Writing them here gives them the one thing they
    # were missing: a committed artefact that says what they are.
    out = ROOT / "results" / "tables" / "release_facts.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("check,value,n,note\n" + "".join(
        f"{k},{v},,derived from the built release by scripts/release_consistency.py\n"
        for k, v in facts))

    if unstated:
        log("")
        log("  NOT STATED ANYWHERE, so the pattern may have stopped matching rather than the")
        log("  documents having stopped claiming it: " + ", ".join(unstated))

    log("")
    if problems:
        log(f"  {len(problems)} STALE RELEASE CLAIM(S):")
        for p in problems:
            log(f"    {p}")
        return 1
    log("  every derived count matches every document that states it")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
