"""D1: where the manuscript repeats itself, measured rather than remembered.

    python scripts/prose_audit.py
    python scripts/prose_audit.py --top 40

A 35-40% length cut needs a target list, and reading for redundancy does not scale over 20,000
words across six sections written at different times. This finds the four kinds of slack that
a cut should take first, in the order they are worth taking:

  1. NUMBERS QUOTED IN MORE THAN ONE SECTION. The strongest signal in this manuscript. A figure
     stated in Results and restated in the Introduction and again in the Discussion is two
     restatements, and each one is also a place the number can drift out of agreement with the
     table. Ranked by how many sections repeat it.
  2. NEAR-DUPLICATE SENTENCES ACROSS SECTIONS, by token overlap. Catches the Introduction
     paraphrasing a Results sentence, which no exact-match scan sees.
  3. THE LONGEST PARAGRAPHS, because a cut works on paragraphs and the largest are where the
     recoverable words are.
  4. SENTENCES OVER 45 WORDS, which are usually two sentences and a subordinate clause.

It reports and does not edit. Every candidate needs a judgement about whether the repetition is
load-bearing -- an abstract is supposed to restate -- and the point is to make that judgement on
a list rather than on a memory of having read the file.
"""

import argparse
import re
from collections import defaultdict
from difflib import SequenceMatcher
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SECT = ROOT / "manuscript" / "sections"
SKIP = {"bibliography.tex"}
STOP = set("the a an and or of to in for by is are was were be been with that this these those "
           "it its as at on from not no than then so such which what where when how we our "
           "us they their there here but if while also both each one two more most less".split())
MIN_SENT_WORDS = 8
LONG_SENT = 45
DUP_RATIO = 0.62


def strip_tex(s):
    """Prose only: drop math, tables, figures, comments and macro names."""
    s = re.sub(r"(?m)^\s*%.*$", " ", s)
    s = re.sub(r"\\begin\{(table|figure|tabular|equation|itemize|enumerate)\*?\}.*?"
               r"\\end\{\1\*?\}", " ", s, flags=re.S)
    s = re.sub(r"\$[^$]*\$", " NUM ", s)
    s = re.sub(r"\\(cite[a-z]*|ref|label|doi|texttt|emph|textbf)\{[^}]*\}", " ", s)
    s = re.sub(r"\\[a-zA-Z]+\*?", " ", s)
    s = re.sub(r"[{}~\\]", " ", s)
    return re.sub(r"\s+", " ", s)


def sentences(text):
    for raw in re.split(r"(?<=[.!?])\s+(?=[A-Z(])", text):
        w = raw.split()
        if len(w) >= MIN_SENT_WORDS:
            yield raw.strip()


def paragraphs(src):
    return [p for p in re.split(r"\n\s*\n", src) if p.strip()]


def tokens(s):
    return {w for w in re.findall(r"[a-z]{3,}", s.lower()) if w not in STOP}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--top", type=int, default=20)
    a = p.parse_args()

    files = sorted(f for f in SECT.glob("*.tex") if f.name not in SKIP)
    prose, sents, paras = {}, {}, {}
    for f in files:
        src = f.read_text()
        prose[f.name] = strip_tex(src)
        sents[f.name] = list(sentences(prose[f.name]))
        paras[f.name] = [(len(strip_tex(q).split()), strip_tex(q)[:90])
                         for q in paragraphs(src)]

    total = sum(len(v.split()) for v in prose.values())
    print(f"=== prose audit: {total} words over {len(files)} sections "
          f"(bibliography excluded) ===\n")
    for f in files:
        n = len(prose[f.name].split())
        print(f"  {f.name:26s} {n:6d} words  {100 * n / total:5.1f}%  "
              f"{len(sents[f.name]):4d} sentences")

    # 1. NUMBERS RESTATED ACROSS SECTIONS.
    where = defaultdict(set)
    for f in files:
        for m in re.findall(r"\d+\.\d{2,}", Path(SECT / f.name).read_text()):
            where[m].add(f.name.replace(".tex", ""))
    multi = sorted(((len(v), k, sorted(v)) for k, v in where.items() if len(v) > 1),
                   reverse=True)
    print(f"\n  {len(multi)} numeric values appear in more than one section. "
          f"Top {a.top}:")
    for n, val, secs in multi[:a.top]:
        print(f"    {val:>12s}  in {n}: {', '.join(secs)}")

    # 2. NEAR-DUPLICATE SENTENCES ACROSS SECTIONS.
    print(f"\n  near-duplicate sentence pairs across sections (token overlap "
          f">= {DUP_RATIO}):")
    pairs = []
    for i, f1 in enumerate(files):
        for f2 in files[i + 1:]:
            for s1 in sents[f1.name]:
                t1 = tokens(s1)
                if len(t1) < 6:
                    continue
                for s2 in sents[f2.name]:
                    t2 = tokens(s2)
                    if len(t2) < 6:
                        continue
                    j = len(t1 & t2) / len(t1 | t2)
                    if j >= DUP_RATIO:
                        pairs.append((j, f1.name, f2.name, s1, s2))
    pairs.sort(reverse=True)
    print(f"    {len(pairs)} pairs")
    for j, f1, f2, s1, s2 in pairs[:a.top]:
        print(f"    [{j:.2f}] {f1} <-> {f2}")
        print(f"        A: {s1[:150]}")
        print(f"        B: {s2[:150]}")

    # 3. LONGEST PARAGRAPHS.
    allp = sorted(((n, f.name, head) for f in files for n, head in paras[f.name]),
                  reverse=True)
    print("\n  longest paragraphs (a cut works on these first):")
    for n, f, head in allp[:a.top]:
        print(f"    {n:4d}w  {f:24s} {head}")

    # 4. OVERLONG SENTENCES.
    longs = sorted(((len(s.split()), f.name, s) for f in files for s in sents[f.name]),
                   reverse=True)
    over = [x for x in longs if x[0] >= LONG_SENT]
    print(f"\n  {len(over)} sentences of {LONG_SENT}+ words:")
    for n, f, s in over[:a.top]:
        print(f"    {n:3d}w  {f:24s} {s[:130]}")

    # SEQUENCE-LEVEL DUPLICATION, as a second opinion on (2). Token overlap misses word order,
    # so a sentence and its own reordering score identically; difflib does not. Reported only
    # as a count, because the pair list above is the actionable form.
    n_seq = 0
    for i, f1 in enumerate(files):
        for f2 in files[i + 1:]:
            for s1 in sents[f1.name]:
                for s2 in sents[f2.name]:
                    if abs(len(s1) - len(s2)) < 40 and \
                            SequenceMatcher(None, s1, s2).ratio() >= 0.72:
                        n_seq += 1
    print(f"\n  {n_seq} cross-section pairs also exceed 0.72 character-level similarity")


if __name__ == "__main__":
    main()
