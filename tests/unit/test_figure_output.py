"""Every shipped figure must be vector, with embedded searchable fonts and legible sizes.

WHY A TEST. These are properties a journal checks and a reader notices, and they are all
silently reversible: one `imshow`, one rasterized=True, one matplotlib default restored, and a
figure that was vector becomes a bitmap or its text becomes unselectable. None of that would
fail any existing gate, and none of it is visible in a screenshot.

TYPE 3 IS THE ONE THAT MATTERS MOST. It is matplotlib's default for PDF, and it embeds glyphs as
PostScript drawing programs rather than as a font: the text is not selectable, not searchable
and not reliably extractable, which also breaks the pypdf reading this project uses to check the
manuscript. `pdf.fonttype: 42` is set in scripts/figures.py for that reason and this asserts it
stays set.
"""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
FIGDIR = ROOT / "manuscript" / "figures"
MIN_FONT_PT = 6.5
MIN_SAVE_DPI = 300


def _figures():
    return sorted(FIGDIR.glob("*.pdf"))


def test_the_submission_package_carries_figures():
    assert _figures(), f"no figure PDFs under {FIGDIR}; run manuscript/build.sh"


@pytest.mark.parametrize("path", _figures(), ids=lambda p: p.name)
def test_figure_is_vector_with_embedded_fonts(path):
    pypdf = pytest.importorskip("pypdf")
    page = pypdf.PdfReader(str(path)).pages[0]
    res = page.get("/Resources", {})

    xobj = res.get("/XObject", {})
    rasters = [k for k in xobj
               if str(xobj[k].get_object().get("/Subtype")) == "/Image"] if xobj else []
    assert not rasters, (
        f"{path.name} embeds {len(rasters)} raster image(s); a figure built from plotted data "
        "should be pure vector. Check for imshow or rasterized=True.")

    fonts = res.get("/Font", {})
    assert fonts, f"{path.name} embeds no font; its text may have been drawn as paths"
    subtypes = {str(fonts[k].get_object().get("/Subtype")) for k in fonts}
    assert "/Type3" not in subtypes, (
        f"{path.name} uses Type 3 fonts, matplotlib's PDF default. Its text is neither "
        "selectable nor searchable and pypdf cannot extract it. Set pdf.fonttype = 42.")


def test_figure_settings_are_print_ready():
    """The two settings the output above depends on, read from the source."""
    src = (ROOT / "scripts" / "figures.py").read_text()

    m = re.search(r"SAVE_DPI\s*=\s*(\d+)", src)
    assert m, "SAVE_DPI not found in scripts/figures.py"
    assert int(m.group(1)) >= MIN_SAVE_DPI, (
        f"SAVE_DPI is {m.group(1)}; journals ask for at least {MIN_SAVE_DPI}")

    assert re.search(r'"pdf\.fonttype":\s*42', src), (
        "pdf.fonttype is not 42; see the Type 3 note above")
    assert re.search(r'"ps\.fonttype":\s*42', src), "ps.fonttype is not 42"


def test_no_font_is_too_small_to_read():
    """A 5 pt tick label survives on screen and disappears in print."""
    src = (ROOT / "scripts" / "figures.py").read_text()
    sizes = [float(x) for x in re.findall(r"fontsize=([0-9.]+)", src)]
    assert sizes, "no explicit fontsize= found; this test would be vacuous"
    assert min(sizes) >= MIN_FONT_PT, (
        f"smallest explicit fontsize is {min(sizes)} pt, below the {MIN_FONT_PT} pt floor")
    m = re.search(r'"font\.size":\s*([0-9.]+)', src)
    assert m and float(m.group(1)) >= MIN_FONT_PT, "the default font.size is below the floor"
