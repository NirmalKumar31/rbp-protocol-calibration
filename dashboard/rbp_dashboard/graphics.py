"""Inline SVG diagrams that explain the method, drawn as page content rather than background.

Every diagram here is a SCHEMATIC. None is plotted from data, none is to scale, and none is a
structural claim about a molecule. They exist because the study's subject, how the unbound
comparison windows were chosen, is a procedure rather than a number, and a reader outside the
field cannot get it from a bar chart. Each one carries a visible "schematic" label for that
reason.

They are content and not decoration, so they are not `aria-hidden`: each takes a `title` that
becomes the accessible name.

Animation is SMIL rather than CSS, so a diagram animates wherever it is embedded without
depending on a stylesheet reaching it. Streamlit injects markdown into containers whose
ancestors carry their own animations, and an ancestor animation creates a containing block that
breaks fixed positioning; keeping motion inside the SVG avoids the whole class of problem.
"""

from __future__ import annotations

import math

from . import theme

STEEL, BRASS, SLATE = theme.STEEL, theme.BRASS, theme.SLATE
INK, MUTED, FAINT = theme.INK, theme.INK_MUTED, theme.INK_FAINT
BASE = theme.NUCLEOTIDE


def _frame(body: str, width: int, height: int, title: str, note: str = "schematic") -> str:
    """Wrap a diagram with its accessible name and its schematic label."""
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'width="100%" height="auto" role="img" style="display:block;max-width:100%">'
        f"<title>{title}</title>"
        f'<text x="{width - 6}" y="13" text-anchor="end" fill="{FAINT}" '
        f'font-family="ui-monospace,monospace" font-size="9" '
        f'letter-spacing="1.4">{note.upper()}</text>'
        f"{body}</svg>"
    )


def hero() -> str:
    """An RNA strand with a protein docked on it, breathing. The subject of the study, once.

    Two sine strands and a lozenge sitting over one stretch of them. The lozenge is the
    RNA-binding protein and the stretch under it is a bound window; everything else along the
    strand is a candidate for the unbound set, which is the whole question.
    """
    w, h, mid = 900, 190, 112
    pts_a, pts_b, rungs = [], [], []
    x = 0.0
    while x <= w:
        t = (x / w) * 4.4 * math.pi
        pts_a.append(f"{x:.1f},{mid + 26 * math.sin(t):.1f}")
        pts_b.append(f"{x:.1f},{mid + 26 * math.sin(t + math.pi):.1f}")
        x += 4
    for i in range(64):
        bx = w * i / 63
        t = (bx / w) * 4.4 * math.pi
        y1, y2 = mid + 26 * math.sin(t), mid + 26 * math.sin(t + math.pi)
        b = "ACGU"[i % 4]
        op = 0.14 + 0.5 * abs(y1 - y2) / 52
        rungs.append(
            f'<line x1="{bx:.1f}" y1="{y1:.1f}" x2="{bx:.1f}" y2="{y2:.1f}" '
            f'stroke="{BASE[b]}" stroke-width="1.5" opacity="{op:.2f}"/>')

    body = (
        f'<g>{"".join(rungs)}</g>'
        f'<polyline points="{" ".join(pts_a)}" fill="none" stroke="{STEEL}" '
        f'stroke-width="2" opacity="0.85" stroke-linecap="round"/>'
        f'<polyline points="{" ".join(pts_b)}" fill="none" stroke="{SLATE}" '
        f'stroke-width="2" opacity="0.7" stroke-linecap="round"/>'
        # the protein, drifting along the strand it binds
        f'<g opacity="0.95">'
        f'<animateTransform attributeName="transform" type="translate" '
        f'values="0 0; 96 0; 0 0" dur="17s" repeatCount="indefinite" '
        f'calcMode="spline" keySplines="0.4 0 0.6 1; 0.4 0 0.6 1" keyTimes="0;0.5;1"/>'
        f'<ellipse cx="300" cy="{mid}" rx="62" ry="30" fill="{BRASS}" opacity="0.13"/>'
        f'<ellipse cx="300" cy="{mid}" rx="62" ry="30" fill="none" stroke="{BRASS}" '
        f'stroke-width="1.6" opacity="0.85">'
        f'<animate attributeName="ry" values="30;33;30" dur="4.5s" repeatCount="indefinite"/>'
        f"</ellipse>"
        f'<text x="300" y="{mid + 4}" text-anchor="middle" fill="{BRASS}" '
        f'font-family="ui-monospace,monospace" font-size="11" letter-spacing="1.2">RBP</text>'
        f"</g>"
        f'<text x="300" y="{mid + 58}" text-anchor="middle" fill="{MUTED}" '
        f'font-family="ui-sans-serif,sans-serif" font-size="11">bound window</text>'
        f'<text x="700" y="{mid + 58}" text-anchor="middle" fill="{FAINT}" '
        f'font-family="ui-sans-serif,sans-serif" font-size="11">'
        f"everything else is a candidate for &quot;unbound&quot;</text>"
        f'<line x1="560" y1="{mid + 44}" x2="840" y2="{mid + 44}" stroke="{FAINT}" '
        f'stroke-width="1" stroke-dasharray="3 4" opacity="0.6"/>'
    )
    return _frame(body, w, h, "An RNA-binding protein bound to one window of an RNA strand")


def negative_sets() -> str:
    """The three protocols, side by side, as the selection procedures they are.

    This is the diagram the dashboard most needed. The result is that three defensible ways of
    filling the bottom row change the measured answer almost fivefold, and that is impossible to
    feel from a bar chart alone.
    """
    w, h = 900, 300
    bound = [(96, 40), (248, 34), (392, 46), (566, 38), (726, 42)]
    rows = [
        ("GC-matched", STEEL, [(150, 40), (300, 34), (440, 46), (612, 38), (772, 42)],
         "pick stretches with a similar G+C fraction"),
        ("Dinucleotide-matched", BRASS, [(140, 40), (292, 34), (430, 46), (604, 38), (764, 42)],
         "also match adjacent letter pairs, a stricter test"),
        ("Bias-aware", SLATE, [(176, 40), (326, 34), (462, 46), (636, 38), (796, 42)],
         "control the assay's own technical bias instead"),
    ]
    body = [
        f'<text x="0" y="30" fill="{MUTED}" font-family="ui-sans-serif,sans-serif" '
        f'font-size="11.5">The transcript, with the windows the protein actually binds</text>',
        f'<line x1="0" y1="52" x2="{w}" y2="52" stroke="{theme.BORDER}" stroke-width="2"/>',
    ]
    for i, (bx, bw) in enumerate(bound):
        body.append(
            f'<rect x="{bx}" y="43" width="{bw}" height="18" rx="2" fill="{BASE["A"]}" '
            f'opacity="0.75"><animate attributeName="opacity" values="0.55;0.85;0.55" '
            f'dur="5s" begin="{i * 0.4}s" repeatCount="indefinite"/></rect>')
    body.append(
        f'<text x="{w}" y="38" text-anchor="end" fill="{BASE["A"]}" '
        f'font-family="ui-monospace,monospace" font-size="9" letter-spacing="1">BOUND</text>')

    for j, (name, colour, picks, caption) in enumerate(rows):
        y = 104 + j * 64
        body.append(
            f'<text x="0" y="{y - 8}" fill="{colour}" '
            f'font-family="ui-monospace,monospace" font-size="9.5" letter-spacing="1.1">'
            f"{name.upper()}</text>")
        body.append(
            f'<text x="150" y="{y - 8}" fill="{FAINT}" '
            f'font-family="ui-sans-serif,sans-serif" font-size="10.5">{caption}</text>')
        body.append(
            f'<line x1="0" y1="{y + 9}" x2="{w}" y2="{y + 9}" stroke="{theme.BORDER_SOFT}" '
            f'stroke-width="1.5"/>')
        for k, (px, pw) in enumerate(picks):
            body.append(
                f'<rect x="{px}" y="{y + 1}" width="{pw}" height="16" rx="2" '
                f'fill="{colour}" opacity="0.5">'
                f'<animate attributeName="opacity" values="0;0.62;0.5" dur="1.1s" '
                f'begin="{j * 0.5 + k * 0.09}s" fill="freeze"/></rect>')
    return _frame("".join(body), w, h,
                  "Three different procedures for choosing unbound comparison windows")


def cross_fitting() -> str:
    """Why the conventional estimator reports signal that cannot exist.

    Five folds. In the two-stage route the score being evaluated was produced by a model that
    saw the evaluation fold; in the cross-fitted route it was not. The diagram is the argument.
    """
    w, h = 900, 250
    fold_w, gap, x0 = 148, 14, 96
    body = []
    for row, (label, colour, leak) in enumerate((
        ("Two-stage", BRASS, True), ("Cross-fitted", STEEL, False),
    )):
        y = 62 + row * 108
        body.append(
            f'<text x="0" y="{y + 20}" fill="{colour}" '
            f'font-family="ui-monospace,monospace" font-size="9.5" letter-spacing="1.1">'
            f"{label.upper()}</text>")
        for i in range(5):
            x = x0 + i * (fold_w + gap)
            held = i == 2
            body.append(
                f'<rect x="{x}" y="{y}" width="{fold_w}" height="34" rx="3" '
                f'fill="{colour if held else theme.RAISED}" '
                f'opacity="{0.28 if held else 1}" stroke="{colour if held else theme.BORDER}" '
                f'stroke-width="1.3"/>')
            body.append(
                f'<text x="{x + fold_w / 2}" y="{y + 22}" text-anchor="middle" '
                f'fill="{colour if held else FAINT}" '
                f'font-family="ui-monospace,monospace" font-size="9.5">'
                f'{"HELD OUT" if held else "TRAIN"}</text>')
        if leak:
            sx = x0 + 2 * (fold_w + gap) + fold_w / 2
            body.append(
                f'<path d="M {x0 + fold_w / 2} {y + 44} Q {sx} {y + 76} {sx} {y + 46}" '
                f'fill="none" stroke="{BRASS}" stroke-width="1.5" stroke-dasharray="4 3" '
                f'opacity="0.9"><animate attributeName="stroke-dashoffset" from="0" to="-14" '
                f'dur="1.1s" repeatCount="indefinite"/></path>')
            body.append(
                f'<text x="{sx + 14}" y="{y + 74}" fill="{BRASS}" '
                f'font-family="ui-sans-serif,sans-serif" font-size="10.5">'
                f"the score was made by a model that saw this fold</text>")
        else:
            body.append(
                f'<text x="{x0}" y="{y + 62}" fill="{MUTED}" '
                f'font-family="ui-sans-serif,sans-serif" font-size="10.5">'
                f"the held-out fold is scored only by models that never saw it</text>")
    body.insert(0,
                f'<text x="0" y="28" fill="{MUTED}" font-family="ui-sans-serif,sans-serif" '
                f'font-size="11.5">The same data, split five ways, scored two ways</text>')
    return _frame("".join(body), w, h, "How cross-fitting closes the outer-fold information path")


def composition() -> str:
    """What the baseline actually is: counting, not understanding.

    The comparison the whole paper rests on is against a model that can only count letters, so
    it is worth showing how little that is.
    """
    w, h = 900, 150
    body = [
        f'<text x="0" y="26" fill="{MUTED}" font-family="ui-sans-serif,sans-serif" '
        f'font-size="11.5">The baseline sees only these 19 numbers. No order, no motifs.</text>'
    ]
    seq = "GAUUACAGAUUACAGGCUAACGUAGCUAAUCGGAUCCAUGGCAUUAGC"
    for i, b in enumerate(seq):
        x = i * 15
        body.append(
            f'<text x="{x}" y="66" fill="{BASE[b]}" font-family="ui-monospace,monospace" '
            f'font-size="13" opacity="0.9">{b}</text>')
    labels = [("A 24%", BASE["A"]), ("C 21%", BASE["C"]), ("G 27%", BASE["G"]),
              ("U 28%", BASE["U"]), ("GC 48%", MUTED), ("+ 14 dinucleotide pairs", FAINT)]
    x = 0
    for i, (text, colour) in enumerate(labels):
        body.append(
            f'<rect x="{x}" y="92" width="{len(text) * 7.6 + 18:.0f}" height="26" rx="3" '
            f'fill="{theme.RAISED}" stroke="{theme.BORDER}" stroke-width="1">'
            f'<animate attributeName="opacity" values="0;1" dur="0.5s" '
            f'begin="{i * 0.12}s" fill="freeze"/></rect>')
        body.append(
            f'<text x="{x + 9}" y="109" fill="{colour}" '
            f'font-family="ui-monospace,monospace" font-size="10.5">'
            f'<animate attributeName="opacity" values="0;1" dur="0.5s" '
            f'begin="{i * 0.12}s" fill="freeze"/>{text}</text>')
        x += len(text) * 7.6 + 32
    return _frame("".join(body), w, h, "The 19 composition features the baseline model uses")


def matching_procedure() -> str:
    """The matching loop, animated: measure the bound window, search the pool, keep the match.

    This is the procedure the whole study varies, and it is a loop rather than a number, so it
    needs a diagram. Four stages play on a repeating cycle so a reader can watch a single pair
    being formed rather than read four captions.
    """
    w, h = 900, 290
    dur = 9.0
    stages = [
        ("1", "take a window the protein binds", 0.0),
        ("2", "measure its letter composition", dur * 0.25),
        ("3", "search the pool of unbound windows", dur * 0.5),
        ("4", "keep it as a negative", dur * 0.75),
    ]
    body = []

    # stage labels, lighting in turn
    box_w, gap = 214, 14
    for i, (num, text, begin) in enumerate(stages):
        x = i * (box_w + gap)
        body.append(
            f'<g><rect x="{x}" y="18" width="{box_w}" height="30" rx="3" fill="{theme.RAISED}" '
            f'stroke="{theme.BORDER}" stroke-width="1"/>'
            f'<rect x="{x}" y="18" width="{box_w}" height="30" rx="3" fill="{STEEL}" opacity="0">'
            f'<animate attributeName="opacity" values="0;0.16;0.16;0" '
            f'dur="{dur}s" begin="{begin}s" keyTimes="0;0.06;0.22;0.3" '
            f'repeatCount="indefinite"/></rect>'
            f'<text x="{x + 12}" y="38" fill="{STEEL}" font-family="ui-monospace,monospace" '
            f'font-size="11">{num}</text>'
            f'<text x="{x + 30}" y="38" fill="{MUTED}" font-family="ui-sans-serif,sans-serif" '
            f'font-size="10.5">{text}</text></g>')

    # the bound window
    body.append(f'<rect x="40" y="96" width="120" height="30" rx="3" fill="{BASE["A"]}" '
                f'opacity="0.68"/>')
    body.append(f'<text x="100" y="116" text-anchor="middle" fill="{INK}" '
                f'font-family="ui-monospace,monospace" font-size="11">BOUND</text>')
    body.append(f'<text x="100" y="146" text-anchor="middle" fill="{FAINT}" '
                f'font-family="ui-monospace,monospace" font-size="10">GC 0.58</text>')

    # the candidate pool, one row of boxes with differing GC
    gcs = [0.31, 0.44, 0.57, 0.72, 0.39, 0.61, 0.48, 0.66]
    for i, gc in enumerate(gcs):
        x = 330 + i * 68
        hit = abs(gc - 0.58) < 0.02
        body.append(
            f'<rect x="{x}" y="96" width="56" height="30" rx="3" '
            f'fill="{theme.BORDER}" opacity="0.55"/>'
            f'<text x="{x + 28}" y="116" text-anchor="middle" fill="{FAINT}" '
            f'font-family="ui-monospace,monospace" font-size="9.5">{gc:.2f}</text>')
        if hit:
            body.append(
                f'<rect x="{x}" y="96" width="56" height="30" rx="3" fill="none" '
                f'stroke="{STEEL}" stroke-width="2" opacity="0">'
                f'<animate attributeName="opacity" values="0;0;1;1;0" dur="{dur}s" '
                f'keyTimes="0;0.5;0.62;0.95;1" repeatCount="indefinite"/></rect>')
            body.append(
                f'<text x="{x + 28}" y="146" text-anchor="middle" fill="{STEEL}" '
                f'font-family="ui-monospace,monospace" font-size="9.5" opacity="0">KEPT'
                f'<animate attributeName="opacity" values="0;0;1;1;0" dur="{dur}s" '
                f'keyTimes="0;0.62;0.72;0.95;1" repeatCount="indefinite"/></text>')

    # the scanning head
    body.append(
        f'<rect x="326" y="92" width="64" height="38" rx="4" fill="none" stroke="{BRASS}" '
        f'stroke-width="1.6" opacity="0.9">'
        f'<animate attributeName="x" values="326;326;530;530" dur="{dur}s" '
        f'keyTimes="0;0.5;0.62;1" repeatCount="indefinite"/>'
        f'<animate attributeName="opacity" values="0;0;0.9;0.9;0" dur="{dur}s" '
        f'keyTimes="0;0.48;0.52;0.6;0.66" repeatCount="indefinite"/></rect>')
    body.append(f'<text x="330" y="88" fill="{FAINT}" font-family="ui-monospace,monospace" '
                f'font-size="9">CANDIDATE POOL, GC FRACTION</text>')

    # the three protocols diverge from the same step
    # Row labels use INK_MUTED, not the series colour. Slate at 9.5px on near-black measures
    # 3.26:1, which is acceptable for a bar and not for type. The swatch beside each label
    # carries the identity instead.
    notes = [("GC-matched", "accept a GC gap under 0.05", STEEL, 0.948),
             ("Dinucleotide-matched", "also match adjacent pairs", BRASS, 0.850),
             ("Bias-aware", "do not match composition at all", SLATE, 0.204)]
    for i, (name, rule, colour, frac) in enumerate(notes):
        y = 196 + i * 30
        body.append(
            f'<rect x="40" y="{y - 8}" width="8" height="8" rx="1.5" fill="{colour}"/>')
        body.append(
            f'<text x="54" y="{y}" fill="{MUTED}" font-family="ui-monospace,monospace" '
            f'font-size="9.5" letter-spacing="0.8">{name.upper()}</text>')
        body.append(
            f'<text x="250" y="{y}" fill="{FAINT}" font-family="ui-sans-serif,sans-serif" '
            f'font-size="10.5">{rule}</text>')
        body.append(
            f'<rect x="560" y="{y - 11}" width="280" height="13" rx="2" '
            f'fill="{theme.BORDER}" opacity="0.4"/>')
        body.append(
            f'<rect x="560" y="{y - 11}" width="0" height="13" rx="2" fill="{colour}" '
            f'opacity="0.7"><animate attributeName="width" values="0;{280 * frac:.0f}" '
            f'dur="1.3s" begin="{0.5 + i * 0.28}s" fill="freeze"/></rect>')
        body.append(
            f'<text x="850" y="{y}" fill="{colour}" font-family="ui-monospace,monospace" '
            f'font-size="9.5">{frac:.0%}</text>')
    body.append(
        f'<text x="560" y="182" fill="{FAINT}" font-family="ui-monospace,monospace" '
        f'font-size="9">SHARE OF PAIRS MATCHED WITHIN 0.05 GC</text>')

    return _frame("".join(body), w, h,
                  "How a bound window is paired with an unbound one, and how the three "
                  "protocols differ in what they accept")
