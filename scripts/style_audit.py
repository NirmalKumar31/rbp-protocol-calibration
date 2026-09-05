"""Stylometric tells of machine-generated prose, counted rather than sensed.

    python scripts/style_audit.py
    python scripts/style_audit.py --show 8

"Does this read as AI-written?" is usually answered by feel, which makes it unfixable: you
cannot tell whether an edit helped. These are the patterns that actually give LLM prose away,
each counted with its instances printed so a specific sentence can be rewritten and the number
watched to fall.

WHAT IS AND IS NOT A TELL HERE. Long sentences are not a tell in a methods paper; UNIFORM
sentence length is. Explanation is not a tell; explaining the obvious immediately after asserting
it is. A colon is not a tell; a colon in every third sentence is. So the measures are mostly
about VARIANCE and RATE, not about the presence of any construction.

The em-dash count is reported and expected to be zero: this author removes them deliberately, so
a nonzero count means one crept back in, not that the prose is machine-like.
"""

import argparse
import re
import statistics as st
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SECT = ROOT / "manuscript" / "sections"
SKIP = {"bibliography.tex"}

# Signposting adverbs and throat-clearing openers. Each is fine once; a rate above roughly one
# per 400 words is the register of a chatbot rather than a journal.
SIGNPOST = ["Importantly", "Notably", "Crucially", "Interestingly", "Significantly",
            "It is worth noting", "It is important to", "It should be noted",
            "In other words", "That said", "Furthermore", "Moreover", "Additionally",
            "Overall", "Ultimately", "In essence", "Simply put", "Indeed"]
# The single most recognisable LLM rhetorical move.
NOT_X_BUT_Y = re.compile(r"\bis not\s+[^.,;]{2,40}?,?\s*(?:but|it is)\s", re.I)
NOT_ONLY = re.compile(r"\bnot only\b[^.]{0,80}?\bbut\b", re.I)
# "rather" is excluded: in this manuscript it is almost always the contrastive "X rather than
# Y", which is a construction and not a hedge. Counting it put 51 false positives at the top of
# the hedge list and hid the four real ones. It is counted separately below as a verbal tic.
HEDGE = ["arguably", "somewhat", "fairly", "quite", "relatively", "potentially",
         "essentially", "largely", "generally", "typically", "may well", "to some extent"]
# Constructions that are individually fine and collectively a signature. A rate under roughly
# one per 300 words means the reader starts hearing the pattern instead of the argument.
TICS = ["rather than", "which is", "so that", "the fact that", "in order to",
        "it is worth", "what matters", "that is", "and not"]


def strip_tex(s):
    s = re.sub(r"(?m)^\s*%.*$", " ", s)
    s = re.sub(r"\\begin\{(table|figure|tabular|equation|itemize|enumerate)\*?\}.*?"
               r"\\end\{\1\*?\}", " ", s, flags=re.S)
    s = re.sub(r"\$[^$]*\$", " N ", s)
    s = re.sub(r"\\(cite[a-z]*|ref|label|doi|texttt|emph|textbf)\{[^}]*\}", " ", s)
    s = re.sub(r"\\[a-zA-Z]+\*?", " ", s)
    s = re.sub(r"[{}~\\]", " ", s)
    return re.sub(r"\s+", " ", s)


def sentences(t):
    return [x.strip() for x in re.split(r"(?<=[.!?])\s+(?=[A-Z(])", t) if len(x.split()) >= 5]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--show", type=int, default=5)
    a = p.parse_args()

    files = sorted(f for f in SECT.glob("*.tex") if f.name not in SKIP)
    text = {f.name: strip_tex(f.read_text()) for f in files}
    allsents = [(f, s) for f in files for s in sentences(text[f.name])]
    words = sum(len(t.split()) for t in text.values())
    lens = [len(s.split()) for _, s in allsents]

    print(f"=== style audit: {words} words, {len(allsents)} sentences ===\n")

    # 1. CADENCE. Human academic prose varies; generated prose regresses to a comfortable
    # middle. Coefficient of variation below ~0.45 is the flat rhythm readers notice without
    # being able to name.
    cv = st.pstdev(lens) / st.mean(lens)
    print(f"  sentence length  mean {st.mean(lens):.1f}  median {st.median(lens)}  "
          f"sd {st.pstdev(lens):.1f}  CV {cv:.3f}")
    print(f"    short (<12w) {sum(n < 12 for n in lens):4d}   "
          f"long (>40w) {sum(n > 40 for n in lens):4d}")
    if cv < 0.45:
        print("    FLAG: cadence is too uniform; add short declarative sentences")

    # 2. SIGNPOSTING RATE.
    hits = Counter()
    for _f, s in allsents:
        for w in SIGNPOST:
            if re.search(rf"\b{re.escape(w)}\b", s):
                hits[w] += 1
    tot = sum(hits.values())
    per = words / tot if tot else float("inf")
    print(f"\n  signposting adverbs: {tot} instances, one per {per:.0f} words")
    for w, n in hits.most_common(8):
        print(f"    {w:22s} {n}")
    if tot and per < 400:
        print("    FLAG: signposting rate is in chatbot register")

    # 3. THE "not X, but Y" MOVE, which is the most recognisable single construction.
    nx = [(f.name, s) for f, s in allsents if NOT_X_BUT_Y.search(s) or NOT_ONLY.search(s)]
    print(f"\n  'not X but Y' / 'not only ... but' constructions: {len(nx)}")
    for f, s in nx[:a.show]:
        print(f"    {f}: {s[:150]}")

    # 4. REPEATED SENTENCE OPENERS. Three words is enough to catch a habit.
    op = Counter(" ".join(s.split()[:3]) for _, s in allsents)
    rep = [(o, n) for o, n in op.most_common(12) if n >= 4]
    print(f"\n  repeated 3-word openers used 4+ times: {len(rep)}")
    for o, n in rep[:8]:
        print(f"    {n:3d}x  {o}")

    # 5. PUNCTUATION STACKING.
    body = " ".join(text.values())
    for ch, name, lim in ((":", "colon", 120), (";", "semicolon", 120)):
        n = body.count(ch)
        print(f"\n  {name}s: {n}, one per {words / max(n, 1):.0f} words"
              + ("   FLAG: stacking clauses" if n and words / n < lim else ""))
    dash = body.count("—") + body.count("--")
    print(f"  em-dashes: {dash}" + ("   FLAG: one crept back in" if dash else "  (deliberate)"))

    # 6. HEDGE DENSITY.
    hh = Counter()
    for w in HEDGE:
        n = len(re.findall(rf"\b{w}\b", body, re.I))
        if n:
            hh[w] = n
    print(f"\n  hedges: {sum(hh.values())}, one per "
          f"{words / max(sum(hh.values()), 1):.0f} words")
    for w, n in hh.most_common(6):
        print(f"    {w:14s} {n}")

    # 6b. VERBAL TICS. A construction repeated every few hundred words is what a reader
    # registers as "samey" even when no individual sentence is wrong.
    print("\n  repeated constructions (rate under 1 per 300 words is audible):")
    for w in TICS:
        n = len(re.findall(rf"\b{re.escape(w)}\b", body, re.I))
        if n >= 12:
            r = words / n
            print(f"    {w:16s} {n:4d}   one per {r:4.0f} words"
                  + ("   FLAG" if r < 300 else ""))

    # 6c. THE MOVES A HUMAN READER ACTUALLY NOTICES. These came from a blind read by someone
    # outside the field, who reported that the standard chatbot vocabulary was absent but that
    # three argumentative GESTURES fired often enough to read as machine cadence. Counting them
    # is the only way to tell whether an edit pass helped.
    print("\n  argumentative gestures (monotony, not vocabulary, is the tell here):")
    NEGATE = re.compile(
        r"(?:,\s*not\s+(?:a|an|the|that|its|his|her|their|only)\b"      # ", not a setting"
        r"|\bis not\s+[^.,;]{2,50}?[,.]\s*(?:It|it) is\b"                # "is not X. It is Y"
        r"|\bnot\s+[^.,;]{2,50}?\s+but\s+)")
    WORTH = re.compile(r"\b(?:is |are |it is |and a )?worth (?:stating|saying|noting|"
                       r"recording|having|making|reporting)\b", re.I)
    CLEFT = re.compile(r"(?:^|\.\s+)What (?:defeats|we|the|this|it|does|makes|it is)\b")
    for name, rx, limit in (("negate-then-assert", NEGATE, 20),
                            ("'worth stating' meta-move", WORTH, 4),
                            ("cleft opener 'What X is...'", CLEFT, 6)):
        found = [(f.name, s_) for f in files for s_ in sentences(text[f.name])
                 if rx.search(s_)]
        print(f"    {name:28s} {len(found):3d}" + ("   FLAG" if len(found) > limit else ""))
        for fn, s_ in found[:a.show]:
            print(f"        {fn}: {s_[:120]}")
    n_so = len(re.findall(r",\s+so\s", body))
    print(f"    {', so' + ' clause chains':28s} {n_so:3d}   one per {words / max(n_so,1):.0f} "
          f"words" + ("   FLAG" if words / max(n_so, 1) < 200 else ""))

    # 7. PARAGRAPH UNIFORMITY. A paper whose paragraphs are all the same size reads as
    # generated even when every sentence is fine.
    plens = []
    for f in files:
        for para in re.split(r"\n\s*\n", f.read_text()):
            n = len(strip_tex(para).split())
            if n >= 25:
                plens.append(n)
    pcv = st.pstdev(plens) / st.mean(plens)
    print(f"\n  paragraphs: {len(plens)}  mean {st.mean(plens):.0f}w  CV {pcv:.3f}"
          + ("   FLAG: too uniform" if pcv < 0.40 else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
