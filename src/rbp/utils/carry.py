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
released evidence, so it is carried forward exactly as committed, and the fact that this run did
not recompute it is reported to the log.

WHAT THIS IS NOT. It is not a cache that can invent a value. If there is no committed table and
no input, the row is absent and `verify.py`'s presence assertion fails, which is correct: that
is a genuinely unanswerable state and it should be loud.

A CARRIED ROW GOES BACK BYTE-IDENTICAL, and the first version of this module got that wrong.
It appended "not recomputed here; the input is absent" to the note, which seemed like honest
provenance and was in fact two defects. The carried table then differed from the committed one
in exactly the environment where carrying happens, so scripts/cache_idempotence.py failed in
CI; and because the suffix was appended to whatever note it read, a second run would append it
again, growing the note without bound. The whole purpose here is that the table does not
change, so the note cannot either. That the row was not recomputed is reported to the LOG, and
the log is not a released artefact.
"""

import pandas as pd

from .log import log


def carried(out_path, check):
    """Return the committed row for `check` as a dict, or None if there is nothing to carry.

    `out_path` is the summary table the caller is about to overwrite. Read before writing. The
    row is returned UNMODIFIED, for the reason the module docstring gives.
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
    log(f"  carried forward {check!r}: not recomputed here, the input is absent")
    return True
