"""Keep a row a `--from-cache` run cannot recompute, instead of dropping it.

THE DEFECT THIS EXISTS TO CLOSE, which occurred three times and was fixed once. Several
summary scripts emit a row only when an optional input is present: a peak file, a GENCODE
index, an unpacked external deposit. In the environment the release is meant to be checked in,
a clean `git archive` with no window store and no raw data, those inputs are absent, so the
documented `--from-cache` command rewrites the committed table WITHOUT the row.

The consequence is not cosmetic. `verify.py` asserts row presence, so a reader who does exactly
what the documentation says takes a repository from 1083/1083 to four failures, and the failure
looks like a broken reproduction rather than a missing optional input. Running the offline
pipeline in a clean export is the whole reproducibility claim, and it was self-corrupting.

This is the same bug commit f3fab95 fixed for baseline_order.py on 2026-09-01, where a
degrading script "wrote an EMPTY table over committed evidence". That fix was applied to one
script and the class was never swept for; window_centring.py, gene_clustered_cv.py and
external_replication.py all still had it a week later.

WHY CARRY FORWARD RATHER THAN REFUSE. Refusing is right for a script whose whole output needs
the missing input, which is what baseline_order.py does. It is wrong here: these scripts
legitimately recompute everything else from committed evidence, and killing the documented
pipeline over one optional row trades a silent corruption for a loud one. The committed value IS
released evidence. So it is carried forward with its provenance stated in the note, and the note
says the row was not recomputed, because a value that silently changes meaning is what this
whole module is about.

WHAT THIS IS NOT. It is not a cache that can invent a value. If there is no committed table and
no input, the row is absent and `verify.py`'s presence assertion fails, which is correct: that
is a genuinely unanswerable state and it should be loud.
"""

import pandas as pd


def carried(out_path, check, note_suffix="not recomputed here; the input is absent"):
    """Return the committed row for `check` as a dict, or None if there is nothing to carry.

    `out_path` is the summary table the caller is about to overwrite. Read before writing.
    """
    if not out_path.exists():
        return None
    try:
        prev = pd.read_csv(out_path)
    except (OSError, pd.errors.ParserError):
        return None
    if "check" not in prev.columns:
        return None
    hit = prev[prev["check"] == check]
    if hit.empty:
        return None
    row = hit.iloc[0].to_dict()
    old = str(row.get("note") or "").strip()
    row["note"] = f"{old}; {note_suffix}" if old else note_suffix
    return {k: ("" if pd.isna(v) else v) for k, v in row.items()}


def emit(out, out_path, check, value, note="", n="", recomputed=True, **extra):
    """Append `check` to `out`, recomputing it or carrying the committed value forward.

    `recomputed` is the caller's answer to "could I actually measure this here?". When it is
    False the committed row is used if one exists, and nothing is appended if it does not.
    """
    if recomputed:
        out.append({"check": check, "value": value, "ci_low": "", "ci_high": "", "n": n,
                    "note": note, **extra})
        return True
    row = carried(out_path, check)
    if row is None:
        return False
    out.append(row)
    return True
