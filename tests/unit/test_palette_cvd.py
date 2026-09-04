"""Every pair of colours that shares a panel must be separable, including under CVD.

WHY THIS IS A TEST AND NOT A JUDGEMENT. The separation is computable, so it is computed. The
first run of this check failed on three pairs, and all three failed at NORMAL vision, which
means full-colour readers could not distinguish them either:

    k-mer vs SpliceBERT              14.8   two mid-lightness cool colours
    dinucleotide arm vs Horlacher-1   9.5   blue against purple
    Figure 5b's two greys             8.7   and they are ADJACENT bars

Thresholds, in OKLab Euclidean distance x100:

    normal vision   >= 15   hard floor. Below this the pair is indistinguishable to everyone
                            and no legend or label repairs it, because the reader still has to
                            match a swatch to a mark.
    dichromatic     >= 8    target, under both protanopia and deuteranopia (Vienot 1999
                            simulation on linear sRGB). 6 to 8 passes only where the panel
                            carries a second encoding; below 6 fails.

Only pairs that actually SHARE A PANEL are checked. A palette-wide all-pairs check would
demand nine mutually distinct hues for colours that never appear together, which buys the
reader nothing and costs the distinctness of the pairs that do co-occur.
"""

import re
from itertools import combinations
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
NORMAL_FLOOR = 15.0
CVD_FLOOR = 6.0
CVD_TARGET = 8.0

# Pairs sharing an axis, read off scripts/figures.py panel by panel.
CO_OCCURRING = [
    ("kmer", "cnn", "splicebert"),          # f9a/b, f13: the three model classes
    ("gc", "dinuc", "neg2"),                # f10a/c: the three protocols
    ("gc", "splicebert"),                   # f15b: worse-or-better colouring
    ("composition", "splicebert"),          # f0b: the two cell-line bars
    ("grey_light", "splicebert"),           # f15a/b: raw against normalised
    ("grey_light", "gc"),                   # f15b: raw against the worse-case colour
    ("gc", "dinuc", "theirs", "grey_mid", "grey_light"),   # f14b: five bars
]
# A second encoding is present, so 6-8 is tolerated for these.
SECONDARY_ENCODED = {("grey_mid", "grey_light")}


def palette():
    """The COLOR dict, parsed from the source so the test cannot drift from the figures."""
    src = (ROOT / "scripts" / "figures.py").read_text()
    m = re.search(r"^COLOR = \{(.*?)\}", src, re.S | re.M)
    assert m, "COLOR dict not found in scripts/figures.py"
    return dict(re.findall(r'"(\w+)":\s*"(#[0-9a-fA-F]{6})"', m.group(1)))


def _linear(hex_colour):
    c = np.array([int(hex_colour[i:i + 2], 16) / 255 for i in (1, 3, 5)])
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def _oklab(lin):
    m1 = np.array([[0.4122214708, 0.5363325363, 0.0514459929],
                   [0.2119034982, 0.6806995451, 0.1073969566],
                   [0.0883024619, 0.2817188376, 0.6299787005]])
    m2 = np.array([[0.2104542553, 0.7936177850, -0.0040720468],
                   [1.9779984951, -2.4285922050, 0.4505937099],
                   [0.0259040371, 0.7827717662, -0.8086757660]])
    return m2 @ np.cbrt(m1 @ lin)


# Vienot, Brettel & Mollon (1999), applied in linear sRGB.
_PROT = np.array([[0.0, 1.05118294, -0.05116099], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
_DEUT = np.array([[1.0, 0.0, 0.0], [0.9513092, 0.0, 0.04866992], [0.0, 0.0, 1.0]])


def delta_e(a, b, vision="normal"):
    la, lb = _linear(a), _linear(b)
    if vision != "normal":
        m = _PROT if vision == "protanopia" else _DEUT
        la, lb = np.clip(m @ la, 0, 1), np.clip(m @ lb, 0, 1)
    return float(np.linalg.norm(_oklab(la) - _oklab(lb)) * 100)


def _pairs():
    seen = []
    for group in CO_OCCURRING:
        for a, b in combinations(group, 2):
            if (a, b) not in seen and (b, a) not in seen:
                seen.append((a, b))
    return seen


@pytest.mark.parametrize(("a", "b"), _pairs())
def test_co_occurring_colours_are_separable(a, b):
    pal = palette()
    assert a in pal and b in pal, f"{a} or {b} missing from COLOR"
    ca, cb = pal[a], pal[b]
    normal = delta_e(ca, cb)
    assert normal >= NORMAL_FLOOR, (
        f"{a} {ca} vs {b} {cb}: normal-vision separation {normal:.1f} < {NORMAL_FLOOR}. "
        "Full-colour readers cannot tell these apart; re-step one of them.")
    worst = min(delta_e(ca, cb, "protanopia"), delta_e(ca, cb, "deuteranopia"))
    floor = CVD_FLOOR if (a, b) in SECONDARY_ENCODED or (b, a) in SECONDARY_ENCODED \
        else CVD_TARGET
    assert worst >= floor, (
        f"{a} {ca} vs {b} {cb}: dichromatic separation {worst:.1f} < {floor}")


def test_every_palette_entry_is_used():
    """An unused entry is a colour nobody validated against anything."""
    src = (ROOT / "scripts" / "figures.py").read_text()
    body = src.split("COLOR = {", 1)[1].split("}", 1)[1]
    unused = [k for k in palette() if f'COLOR["{k}"]' not in body]
    assert unused == [], f"unused COLOR entries: {unused}"
