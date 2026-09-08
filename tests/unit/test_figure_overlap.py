"""The overlap detector's own geometry, on synthetic boxes.

Why this exists separately from the checker. `scripts/figure_overlap.py` needs `pdftotext`,
which is not present in every environment, so a suite that only ran it end to end would skip
wherever poppler is absent and the detector's logic would be gated by nothing. These tests
exercise `overlapping()` directly, so the arithmetic is covered everywhere and only the PDF
extraction depends on the external tool.

The case that matters is the real one: `f10_three_protocols` drew panel b's label inside panel
a's title, at 100% of the smaller box, in every build it ever had.
"""

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

fo = pytest.importorskip("figure_overlap")


def box(x0, y0, x1, y1, t):
    return (x0, y0, x1, y1, t)


def test_a_label_drawn_inside_a_title_is_caught():
    """The shipped defect, reduced: a one-character label wholly inside a wider word."""
    ws = [box(0, 0, 40, 10, "5.4x"), box(20, 2, 26, 8, "b")]
    hits = fo.overlapping(ws, 0.30)
    assert len(hits) == 1 and hits[0][2] == pytest.approx(1.0), (
        "a label entirely inside another word must be reported at 100% of the smaller box")


def test_words_that_merely_touch_are_not_a_collision():
    """Adjacent glyph boxes share an edge. Zero-area intersection is not overlap."""
    assert fo.overlapping([box(0, 0, 10, 10, "a"), box(10, 0, 20, 10, "b")], 0.30) == []


def test_a_kerned_sliver_is_below_the_threshold():
    """Two boxes sharing 10% of the smaller one are normal typesetting, not a defect."""
    assert fo.overlapping([box(0, 0, 10, 10, "a"), box(9, 0, 19, 10, "b")], 0.30) == []
    assert len(fo.overlapping([box(0, 0, 10, 10, "a"), box(9, 0, 19, 10, "b")], 0.05)) == 1, \
        "the threshold must be the thing that decides, so a stricter one must see it"


def test_boxes_apart_in_y_do_not_collide_however_close_in_x():
    """Two panel titles at the same x on different rows are the common false positive."""
    assert fo.overlapping([box(0, 0, 50, 10, "title a"), box(0, 20, 50, 30, "title b")],
                          0.30) == []


def test_the_checker_refuses_to_pass_over_zero_figures(tmp_path, monkeypatch):
    """A checker that finds nothing must fail, not report success. This repository has shipped
    two gates that passed having examined an empty list."""
    monkeypatch.setattr(fo, "FIGS", tmp_path)
    monkeypatch.setattr(sys, "argv", ["figure_overlap.py"])
    with pytest.raises(SystemExit) as e:
        fo.main()
    assert "refusing" in str(e.value)
