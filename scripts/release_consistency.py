"""FIXING, and why it is positional.

`--fix` rewrites every stale claim in place. It does NOT search for the number: the scan above
already recorded the exact character offset of each stale value, and the rewrite uses that
offset, right to left within a line so earlier edits cannot shift later ones. It asserts the
bytes at that offset are still the digits it found.

That distinction is the reason this mode exists. Syncing these counts by hand meant editing
about twenty places per change, and doing it with a regex over the digits is what rewrote
fourteen unrelated scientific numbers across the manuscript on 2026-09-06 and 2026-09-07: in a
regex `\\b58\\b` matches the `58` inside `6.58`, because `.` is not a word character. A
positional rewrite cannot make that mistake, because it only ever touches a span the scanner
itself identified as an instance of the fact being synced.

A partial run refuses to fix, for the same reason it refuses to rewrite release_facts.csv.

Counts a release document states about the release must be derived from the release.

    python scripts/release_consistency.py

Why this exists. An external review read the repository and found eight stale figures in the
submission package: 28 pages against a 48-page PDF, 20 references against 26, six main figures
against seven, Tables 1 to 7 against fourteen, 768 numeric assertions against 937, 696 in the
Zenodo template, and a recomputation error of 2.2e-16 against the manuscript's 3.3e-16. Every
one had been correct when it was typed. Every one was copied by hand into a file that no gate
read.

Why the existing gates did not catch them. scripts/audit_manuscript.py scans prose for numbers
with no source, and it now scans these documents too, but it cannot catch these: 28, 20 and 7
are small integers that collide with something in a 855-value haystack, and 2.2e-16 is written
in an exponent form its regex does not match. Orphan-hunting asks "is this number real
somewhere". The question here is the different and stricter one: "is this number still true of
THIS artefact". That needs the artefact, not a haystack.

So the facts are derived and the prose is matched against them. Each fact below is computed
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
# Every public, citable or operational file, not only the prose ones. This listed the five
# markdown documents and the TeX, and a review then found stale "937" counts in three places in
# .github/workflows/ci.yml, one in docs/COST.md and one in a pyproject comment, plus a stale
# span in CITATION.cff -- every one of them outside the scan. The lesson is the same one that
# produced this script: a checker certifies the files it looked at, and the ones it does not
# look at are where the stale numbers go and stay.
#
# A workflow file and a pyproject comment are not prose, but a reader acts on them and a
# citation file is the most quoted artefact in the repository.
DOCS = [ROOT / "README.md",
        ROOT / "docs" / "REPRODUCE.md", ROOT / "docs" / "PANELS.md",
        ROOT / "docs" / "COST.md",
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


def _verify_summary(key):
    """One row of the summary the verifier writes.

    Read from the artefact and not from verify.py's source, because the count is a property
    of a run: gates skip when a table is absent, and a number taken from the source would be
    the number the code COULD reach rather than the one it did.
    """
    f = ROOT / "results" / "tables" / "verify_summary.csv"
    if not f.exists():
        return None
    for ln in f.read_text().strip().splitlines()[1:]:
        k, _, v = ln.partition(",")
        if k.strip() == key:
            return int(v)
    return None


def n_verify_checks():
    return _verify_summary("assertions")


# THE SPLIT, not just the total. README states how many of the 1156 assertions belong to this
# paper rather than to the earlier variant-scoring study, which is the number a reader of THIS
# paper actually wants. Only the total was ever derived, so when a gate was added the total was
# synced and the split was not: the README said 1015 while the verifier printed 1018. A fact
# stated in prose and derived by nothing is the definition of what goes stale here.
def n_paper_assertions():
    return _verify_summary("assertions_this_paper")


def n_legacy_assertions():
    return _verify_summary("assertions_legacy_study")


def n_harness_assertions():
    """The pair that made the README's split fail to add up: 1018 + 136 is 1154, not 1156."""
    return _verify_summary("assertions_harness")


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
    # AN ESCAPED DOLLAR IS NOT A MATHS DELIMITER, and this treated it as one. Writing a cost
    # of "\\$115" in the abstract made the regex below open a maths span at that dollar and close
    # it at the next real one, swallowing 165 words of prose as a single token: the count went
    # from 465 to 300 and the release check reported a stale claim in a document that was not
    # stale. Escaped dollars are taken out of play first.
    body = m.group(1).replace(r"\$", "\x00")
    body = re.sub(r"\$[^$]*\$", " X ", body)            # one token per maths span
    body = body.replace("\x00", " $ ")                  # and a literal dollar is one word
    body = re.sub(r"\\[a-zA-Z]+\*?", " ", body)         # control sequences are not words
    body = re.sub(r"[{}~\\]", " ", body)
    return len([w for w in body.split() if re.search(r"[A-Za-z0-9]", w)])


def collected_tests():
    """Tests pytest actually collects here, not a number typed into the README.

    This count moves when a script is added, not only when a test is. Two suites parametrise
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
    # Summed from the per-file lines, not read off a summary line. pyproject sets addopts="-q"
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
    # Patterns that missed known-stale claims. An audit found "706 tests pass" in the workflow
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
        # "checks" is NOT required after either form. README states the count as a bare
        # comment, `# 1153/1153`, and another wrapped "All 1153 verification" onto the next
        # line, so both escaped --fix and had to be hand-edited on three consecutive rounds.
        # A syncing tool that cannot reach a claim makes that claim the one that goes stale.
        r"All (\d{3,4}) verification",
        r"(\d{3,4})/\d{3,4} checks",
        r"#\s*(\d{3,4})/\d{3,4}\b",
        # BOTH SIDES of an N/N pair. Capturing only the left one rewrote "# 1153/1153" to
        # "# 1156/1153", which is worse than leaving it alone: a self-contradicting line that
        # still looks synced. The offsets differ so the positional fixer handles both.
        r"#\s*\d{3,4}/(\d{3,4})\b",
        r"\d{3,4}/(\d{3,4}) checks",
        r"one command, (\d{3,4}) checks",
        # Hyphenated and singular: "the 1141-assertion verification harness" sat 15 behind
        # and no pattern here could reach it, the same failure as a document nobody scans.
        r"(\d{3,4})-assertion",
    ]),
    "paper assertions": (n_paper_assertions, [
        r"(\d{3,4}) belong to this paper",
    ]),
    "legacy assertions": (n_legacy_assertions, [
        r", (\d{2,4}) to an earlier",
    ]),
    "harness assertions": (n_harness_assertions, [
        r"and (\d+) are\s+the harness checking itself",
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
# This set turned CI red and that was the right bug to have, handled the wrong way. A previous
# audit said an unavailable required fact must fail rather than pass silently, which is correct.
# Making it fail unconditionally was not: the CPU workflow installs docker/requirements-cpu.txt,
# which has no torch by design, so the full suite cannot be COLLECTED there and the job went red
# on an environment limitation rather than on a stale document.
#
# The distinction that actually matters is between "this environment cannot derive the fact" and
# "this environment could derive it but a dependency is missing". The first is a legitimate
# partial run and is reported; the second is a broken environment and fails. --require-all makes
# every fact mandatory, and the full-suite CI job passes it.
REQUIRED = {"tests collected", "verify checks", "paper assertions", "legacy assertions",
            "harness assertions"}


def title_equality():
    """The paper's title, as three artefacts state it. A citation importer reads the CFF.

    Not a count, so it does not fit the FACTS table, and it is exactly the class of drift that
    table cannot catch: for a while the PDF said "depends strongly on" while both CITATION.cff
    titles still said "is set by", so importing the citation produced a different paper.
    """
    import re as _re

    def norm(s):
        return _re.sub(r"\s+", " ", s).strip().rstrip(".").lower()

    tex = (MANUSCRIPT / "paper.tex").read_text()
    m = _re.search(r"\\title\{\\bfseries (.*?)\}\n", tex, _re.S)
    if not m:
        return []
    paper = norm(_re.sub(r"[%\\]", " ", m.group(1)))
    out = []
    cff = ROOT / "CITATION.cff"
    if cff.exists():
        y = __import__("yaml").safe_load(cff.read_text())
        pref = (y.get("preferred-citation") or {}).get("title", "")
        if pref and norm(pref) != paper:
            out.append("CITATION.cff preferred-citation title differs from the manuscript "
                       f"title:\n      cff:   {norm(pref)}\n      paper: {paper}")
    pdf = MANUSCRIPT / "paper.pdf"
    if pdf.exists():
        try:
            import pypdf
            meta = pypdf.PdfReader(str(pdf)).metadata or {}
            got = norm(str(meta.get("/Title", "")))
            if got and got != paper:
                out.append("PDF /Title differs from the manuscript title:\n"
                           f"      pdf:   {got}\n      paper: {paper}")
        except ImportError:
            pass
    return out


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--require-all", action="store_true",
                    help="fail when any fact in REQUIRED cannot be derived here")
    ap.add_argument("--fix", action="store_true",
                    help="rewrite every stale claim in place, positionally; see FIXING below")
    a = ap.parse_args(argv)
    require_all, fixing = a.require_all, a.fix
    if fixing and (a.require_all is False):
        pass                                   # --fix is usable with or without --require-all
    problems, unstated, facts, skipped_required, broken = [], [], [], [], []
    edits = {}                                 # doc -> {line index: [(start, old_len, new)]}
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
                        # Recorded by POSITION, which is the whole point of --fix. The scanner
                        # already knows the exact offset of this number, so rewriting it cannot
                        # touch any other number in the file. Doing the same job with a regex
                        # over the digits is what rewrote fourteen unrelated scientific values
                        # on 2026-09-06 and 09-07: `\b58\b` matches inside `6.58`.
                        edits.setdefault(doc, {}).setdefault(i - 1, []).append(
                            (start, len(str(got)), str(truth)))
        log(f"  {name:22} {truth:<8} stated in {seen} place(s)")
        if not seen:
            unstated.append(name)

    # The facts are committed, which is what makes them quotable. scripts/audit_manuscript.py
    # flags any number in a released document that no table can source, and a page count or a
    # test census lives in no result table -- so widening that audit to the release documents
    # correctly reported them as unsourced. Writing them here gives them the one thing they
    # were missing: a committed artefact that says what they are.
    # A partial run must not rewrite the committed facts. This file is environment-dependent:
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

    # A required fact that could not be derived is a failure, not a pass. `tests collected`
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

    titles = title_equality()
    if titles:
        log("")
        log("  TITLE MISMATCH between artefacts a citation importer reads:")
        for x in titles:
            log(f"    {x}")

    if fixing and edits:
        if skipped_required or broken:
            log("\n  NOT fixing: this environment could not derive every fact, so a rewrite "
                "would encode a partial run's numbers into the documents")
        else:
            n = 0
            for doc, by_line in edits.items():
                lines = doc.read_text().splitlines(keepends=True)
                for idx, spans in by_line.items():
                    # Right to left, so an earlier edit cannot shift a later offset.
                    for start, old_len, new_text in sorted(spans, reverse=True):
                        ln = lines[idx]
                        assert ln[start:start + old_len].isdigit(), (
                            f"{doc}:{idx + 1} offset {start} is not the number it was found at")
                        lines[idx] = ln[:start] + new_text + ln[start + old_len:]
                        n += 1
                doc.write_text("".join(lines))
            log(f"\n  --fix rewrote {n} stale claim(s) in {len(edits)} file(s), by position")
            log("  Re-run without --fix to confirm, and rebuild the PDFs if a .tex changed.")
            problems = []

    log("")
    if problems or broken or titles or (skipped_required and require_all):
        log(f"  {len(problems)} STALE RELEASE CLAIM(S):")
        for p in problems:
            log(f"    {p}")
        return 1
    log("  every derived count matches every document that states it")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
