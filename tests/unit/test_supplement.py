"""The supplement is a document with legends, and the S-numbering cannot drift from the files.

Why this exists. The ten supplementary figures shipped for weeks as loose PDFs named f0 to f8
and f13, with no S-numbers and no legends, while SUBMISSION.md described that state accurately
and an audit correctly said accuracy about an unpublishable supplement is not a supplement.
manuscript/supplementary.tex now carries a legend per figure. That introduces two new ways to
be wrong, and both are mechanical:

  * the S-number to filename mapping is a hand-written table, so it can name a figure that does
    not exist, or omit one that does;
  * a figure can be in results/figures/ and cited by neither paper.tex nor the supplement, which
    is how five figures were orphaned once before, after a section was cut.

So the mapping is parsed out of the LaTeX and checked against the tree, and the partition of
figures into main and supplementary is required to be exact: every built figure is cited exactly
once, by one document or the other.
"""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
MAN = ROOT / "manuscript"
SUP = MAN / "supplementary.tex"
FIGDIR = ROOT / "results" / "figures"


def _mapping():
    """S-number -> filename stem, parsed from the mapping table in supplementary.tex."""
    rows = re.findall(r"Figure S(\d+)\s*&\s*\\texttt\{([^}]+)\}", SUP.read_text())
    return {int(n): f.replace("\\_", "_").removesuffix(".pdf") for n, f in rows}


# THE EXTENSION IS OPTIONAL IN LaTeX AND THE SECTIONS OMIT IT. A first version of this file
# required "\.pdf" and so found zero main figures, which made the partition test report every
# main figure as an orphan. Match with or without it, exactly as build.sh's own glob does.
_INC = re.compile(r"includegraphics\[[^\]]*\]\{figures/([a-z0-9_]+?)(?:\.pdf)?\}")


def _included():
    """Figure stems actually \\includegraphics'd by the supplement."""
    return set(_INC.findall(SUP.read_text()))


def _main_figures():
    text = "\n".join(p.read_text() for p in (MAN / "sections").glob("*.tex"))
    return set(_INC.findall(text))


def test_the_supplement_exists_at_all():
    assert SUP.exists(), "manuscript/supplementary.tex is the supplement; it must be in the tree"


def test_every_mapped_figure_exists_in_the_tree():
    for n, stem in sorted(_mapping().items()):
        assert (FIGDIR / f"{stem}.pdf").exists(), (
            f"the mapping table names Figure S{n} as {stem}.pdf, which is not in results/figures")


def test_the_s_numbers_are_contiguous_from_one():
    keys = sorted(_mapping())
    assert keys == list(range(1, len(keys) + 1)), (
        f"S-numbers must run 1..N with no gaps or repeats; got {keys}")


def test_every_mapped_figure_is_actually_included():
    """A row in the table that no \\includegraphics backs is a promise, not a figure."""
    mapped = set(_mapping().values())
    inc = _included()
    assert mapped == inc, (
        f"mapped but not included: {sorted(mapped - inc)}; "
        f"included but not mapped: {sorted(inc - mapped)}")


def test_main_and_supplementary_partition_the_built_figures():
    """Exact partition. A figure cited by neither document is an orphan, and five were once."""
    built = {p.stem for p in FIGDIR.glob("f*.pdf")}
    main, sup = _main_figures(), _included()
    assert not (main & sup), f"cited by BOTH documents: {sorted(main & sup)}"
    assert main | sup == built, (
        f"cited by neither: {sorted(built - (main | sup))}; "
        f"cited but not built: {sorted((main | sup) - built)}")


@pytest.mark.parametrize("field", ["n}:", "Uncertainty}:", "Source}:"])
def test_every_legend_states_the_reader_facing_essentials(field):
    """Sample size, what the uncertainty is, and which committed table the values come from.

    Checked as a count against the number of figures rather than per caption, because a caption
    is one LaTeX group and splitting them reliably is more machinery than the property needs.
    """
    text = SUP.read_text()
    n_figs = len(_included())
    got = text.count("\\textbf{" + field)
    assert got >= n_figs, (
        f"{got} legends state '{field}' for {n_figs} figures; every figure needs it")


def test_the_built_supplement_carries_every_figure_label():
    """Skips where the PDF is unbuilt, which is the container image and a fresh clone."""
    pdf = MAN / "supplementary.pdf"
    if not pdf.exists():
        pytest.skip("supplementary.pdf not built; run manuscript/build.sh")
    pypdf = pytest.importorskip("pypdf")
    txt = "\n".join(p.extract_text() for p in pypdf.PdfReader(str(pdf)).pages)
    for n in sorted(_mapping()):
        assert f"Figure S{n}" in txt, f"Figure S{n} has no label in the rendered supplement"


def test_the_submission_index_does_not_still_call_the_supplement_unpackaged():
    """The claim that was true and is not any more.

    SUBMISSION.md said the supplementary figures "are loose PDFs and are not yet a packaged
    supplement". Leaving that sentence after building the supplement would be the same defect
    as the fold-leakage disclosure that outlived its defect by three release candidates.
    """
    text = (ROOT / "SUBMISSION.md").read_text()
    assert "not yet a packaged supplement" not in text, (
        "SUBMISSION.md still describes the supplement as unpackaged")
    assert "supplementary.pdf" in text, "SUBMISSION.md must name the built supplement"
