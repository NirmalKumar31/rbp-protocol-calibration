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
# EVERY PUBLIC, CITABLE OR OPERATIONAL FILE, not only the prose ones. This listed the five
# markdown documents and the TeX, and a review then found stale "937" counts in three places in
# .github/workflows/ci.yml, one in docs/COST.md and one in a pyproject comment, plus a stale
# span in CITATION.cff -- every one of them outside the scan. The lesson is the same one that
# produced this script: a checker certifies the files it looked at, and the ones it does not
# look at are where the stale numbers go and stay.
#
# A workflow file and a pyproject comment are not prose, but a reader acts on them and a
# citation file is the most quoted artefact in the repository.
DOCS = [ROOT / "README.md", ROOT / "SUBMISSION.md",
        ROOT / "docs" / "REPRODUCE.md", ROOT / "docs" / "PANELS.md",
        ROOT / "docs" / "ZENODO.md", ROOT / "docs" / "COST.md",
        ROOT / "docs" / "AUDIT-RESPONSE.md",
        ROOT / "CITATION.cff", ROOT / "CHANGELOG.md", ROOT / "pyproject.toml",
        ROOT / ".github" / "workflows" / "ci.yml"] + TEX


def _tex(paths=None):
    return "\n".join(p.read_text() for p in (paths or TEX))


class MissingTool(RuntimeError):
    """The artefact is present but the thing that reads it is not. Not the same as absent."""


def pdf_pages():
    if not PDF.exists():
        return None
    try:
        import pypdf
    except ImportError as e:
        # The PDF IS here; we simply cannot read it. Reporting that as "artefact absent" is how
        # CI silently stopped checking the page count of a file sitting in the checkout.
        raise MissingTool("paper.pdf is present but pypdf is not installed") from e
    return len(pypdf.PdfReader(str(PDF)).pages)


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
    # PATTERNS THAT MISSED KNOWN-STALE CLAIMS. An audit found "706 tests pass" in the workflow
    # header and "716 tests" in the changelog, both invisible to the two narrow forms this
    # started with. A regex set that only matches the phrasings already in the repository
    # certifies the phrasings, not the facts.
    "tests collected": (collected_tests, [
        r"#\s*(\d+),? (?:passed|pass|needs torch)",
        r"tests/\s+(\d+) tests",
        r"\b(\d{3,4}) (?:collected )?tests\b",
        r"\b(\d{3,4}) tests? pass",
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
        r"(\d+) figure environments",
    ]),
    "table environments": (lambda: n_environments("table"), [
        r"Tables 1 to (\d+) are typeset",
        r"(\d+) table environments",
        r"\d+ pages, (\d+) tables",
    ]),
    "verify checks": (n_verify_checks, [
        # Deliberately loose on the adjective. The narrow form missed "937 published
        # assertions" in three files and "937 checks" in a fourth, all of which a reader
        # acts on exactly as if they said "numeric assertions".
        r"(\d{3,4})\s+\w*\s*assertions",
        r"All (\d{3,4}) verification checks",
        r"(\d{3,4})/\d{3,4} checks",
        r"one command, (\d{3,4}) checks",
    ]),
    "abstract words": (abstract_words, [
        r"(\d+) words, no markup",
    ]),
}

# A word count is not a checksum: LaTeX macro stripping and a form's own counter will differ by
# a word or two, and failing a release on that would train everyone to ignore this script.
TOLERANCE = {"abstract words": 4}

# Facts whose absence means the run proved less than it claims.
#
# THIS SET TURNED CI RED AND THAT WAS THE RIGHT BUG TO HAVE, HANDLED THE WRONG WAY. A previous
# audit said an unavailable required fact must fail rather than pass silently, which is correct.
# Making it fail unconditionally was not: the CPU workflow installs docker/requirements-cpu.txt,
# which has no torch by design, so the full suite cannot be COLLECTED there and the job went red
# on an environment limitation rather than on a stale document.
#
# The distinction that actually matters is between "this environment cannot derive the fact" and
# "this environment could derive it but a dependency is missing". The first is a legitimate
# partial run and is reported; the second is a broken environment and fails. --require-all makes
# every fact mandatory, and the full-suite CI job passes it.
REQUIRED = {"tests collected", "verify checks"}


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--require-all", action="store_true",
                    help="fail when any fact in REQUIRED cannot be derived here")
    require_all = ap.parse_args(argv).require_all
    problems, unstated, facts, skipped_required, broken = [], [], [], [], []
    for name, (derive, patterns) in FACTS.items():
        try:
            truth = derive()
        except MissingTool as e:
            log(f"  {name:22} BROKEN ({e})")
            broken.append(f"{name}: {e}")
            continue
        if truth is None:
            log(f"  {name:22} SKIP (cannot be derived in this environment)")
            if name in REQUIRED:
                skipped_required.append(name)
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
    # A PARTIAL RUN MUST NOT REWRITE THE COMMITTED FACTS. This file is environment-dependent:
    # the no-torch CI job cannot derive the test census, so it dropped that row, and the
    # git-diff gate then failed on a table that was correct for the environment that wrote it.
    # The gate was right and the writer was wrong. Only a run that derived everything may write.
    out = ROOT / "results" / "tables" / "release_facts.csv"
    if skipped_required or broken:
        log(f"  not rewriting {out.name}: this environment could not derive "
            f"{', '.join(skipped_required + [b.split(':')[0] for b in broken])}")
    else:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("check,value,n,note\n" + "".join(
            f"{k},{v},,derived from the built release by scripts/release_consistency.py\n"
            for k, v in facts))

    # A REQUIRED FACT THAT COULD NOT BE DERIVED IS A FAILURE, NOT A PASS. `tests collected`
    # returns None wherever torch is absent, which is exactly the CPU environment CI runs in, so
    # the script printed SKIP and exited zero on the one machine whose job is to catch this.
    # Named explicitly rather than "everything must derive", because a manuscript that has not
    # been built genuinely has no page count and failing on that would be noise.
    if broken:
        log("")
        log("  A DEPENDENCY IS MISSING, so a fact that IS derivable here went unchecked:")
        for s in broken:
            log(f"    {s}")

    if skipped_required:
        log("")
        log("  REQUIRED FACTS COULD NOT BE DERIVED, so this run certifies less than it looks:")
        for s in skipped_required:
            log(f"    {s}")
        if require_all:
            log("  --require-all was given, so this is a failure.")
        else:
            log("  Not a failure without --require-all: this environment is a declared subset.")
            log("  The full-suite CI job passes --require-all and is where these are enforced.")

    if unstated:
        log("")
        log("  NOT STATED ANYWHERE, so the pattern may have stopped matching rather than the")
        log("  documents having stopped claiming it: " + ", ".join(unstated))

    log("")
    if problems or broken or (skipped_required and require_all):
        log(f"  {len(problems)} STALE RELEASE CLAIM(S):")
        for p in problems:
            log(f"    {p}")
        return 1
    log("  every derived count matches every document that states it")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
