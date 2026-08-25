"""Dinucleotide-preserving shuffle, for the strongest negative control we can build.

Our matched negatives control region and GC content. The EDA showed that is not enough:
GC constrains G+C but not how those bases are arranged, and TARDBP's positives carry GU
dinucleotides at 3.8x the rate of its negatives. So a model could score well on "G-rich,
C-poor, GU-heavy" without recognising any motif.

A dinucleotide-preserving shuffle destroys motifs while holding mononucleotide AND
dinucleotide frequencies exactly constant. Whatever AUROC survives against these
negatives is genuinely motif-driven rather than compositional.

The method is Altschul and Erikson (1985). A sequence is an Eulerian path through a
multigraph whose vertices are nucleotides and whose edges are the dinucleotides. Any
other Eulerian path through the same multigraph is a sequence with identical dinucleotide
counts. So: shuffle the edge lists, check the arrangement still admits an Eulerian path,
and walk it.

The connectivity condition is the part that is easy to get wrong. After shuffling, the
LAST outgoing edge of every vertex must lead, by repeated following, to the terminal
vertex. If it does not, the walk strands itself with edges unused.
"""

from collections import defaultdict

import numpy as np


def _last_edges_reach(edges, verts, last):
    """Does following each vertex's final edge always arrive at `last`?

    This is the Eulerian-path condition. Without it the walk can dead-end partway and
    silently return a truncated sequence.
    """
    for v in verts:
        if v == last or not edges[v]:
            continue
        seen = {v}
        cur = v
        while cur != last:
            if not edges[cur]:
                return False
            nxt = edges[cur][-1]
            if nxt in seen:            # a cycle that never reaches `last`
                return False
            seen.add(nxt)
            cur = nxt
    return True


def dinuc_shuffle(seq, rng, max_tries=200):
    """A random sequence with exactly the same mono- and dinucleotide counts as `seq`.

    Returns None if no valid arrangement was found, which the caller must handle rather
    than silently accepting a broken sequence. In practice this is rare for 101-nt
    windows over a 4-letter alphabet.
    """
    if len(seq) < 3:
        return seq

    first, last = seq[0], seq[-1]
    edges = defaultdict(list)
    for a, b in zip(seq, seq[1:]):
        edges[a].append(b)
    verts = list(edges)

    for _ in range(max_tries):
        for v in verts:
            rng.shuffle(edges[v])
        if _last_edges_reach(edges, verts, last):
            break
    else:
        return None

    pos = {v: 0 for v in verts}
    out = [first]
    cur = first
    for _ in range(len(seq) - 1):
        if pos[cur] >= len(edges[cur]):
            return None                # should not happen once connectivity holds
        nxt = edges[cur][pos[cur]]
        pos[cur] += 1
        out.append(nxt)
        cur = nxt
    return "".join(out)


def mono_shuffle(seq, rng):
    """Plain shuffle: preserves mononucleotide counts only.

    Kept as a weaker rung on the control ladder. The gap between mono- and
    dinucleotide-shuffled performance is exactly how much a model is exploiting
    dinucleotide structure rather than single-base composition.
    """
    a = np.array(list(seq))
    rng.shuffle(a)
    return "".join(a)


def counts(seq, k=2):
    """k-mer counts, for verifying a shuffle preserved what it claims to."""
    out = defaultdict(int)
    for i in range(len(seq) - k + 1):
        out[seq[i:i + k]] += 1
    return dict(out)


def verify(original, shuffled, k=2):
    """Check the shuffle preserved k-mer counts and changed the arrangement."""
    return {
        "length_same": len(original) == len(shuffled),
        "mono_same": counts(original, 1) == counts(shuffled, 1),
        "dinuc_same": counts(original, 2) == counts(shuffled, 2),
        "differs": original != shuffled,
        "identity": round(sum(a == b for a, b in zip(original, shuffled)) / len(original), 4),
    }


def shuffled_negatives(positives, seed=7, method="dinuc", max_identity=0.9):
    """One shuffled negative per positive sequence.

    Returns (sequences, dropped). A shuffle that comes back almost identical to its
    source is rejected, because a near-copy of a positive is not a negative -- with a
    low-complexity input such as a poly-U tract, a dinucleotide-preserving shuffle can
    legitimately return something very close to the original.
    """
    rng = np.random.default_rng(seed)
    fn = dinuc_shuffle if method == "dinuc" else mono_shuffle
    out, dropped = [], {"failed": 0, "too_similar": 0}
    for s in positives:
        best = None
        for _ in range(10):
            cand = fn(s, rng) if method == "dinuc" else fn(s, rng)
            if cand is None:
                continue
            ident = sum(a == b for a, b in zip(s, cand)) / len(s)
            if ident <= max_identity:
                best = cand
                break
            best = best or cand
        if best is None:
            dropped["failed"] += 1
            out.append(None)
        elif sum(a == b for a, b in zip(s, best)) / len(s) > max_identity:
            dropped["too_similar"] += 1
            out.append(None)
        else:
            out.append(best)
    return out, dropped
