"""Line-buffered progress output, in one place.

Fifty scripts each carried an identical three-line `log`. That is harmless until someone wants
to change how progress is reported -- add a timestamp, route to a file, silence it under a
flag -- at which point it is fifty edits, and the copies drift instead. Living in `rbp.utils`
rather than in a `scripts/` helper means it is importable from anywhere `src` is on the path,
which every entry point already arranges.

flush=True is the point of the function. These scripts are run under `nohup` and read through a
log file while they work, and Python's default block buffering makes a long run look hung for
minutes at a time.
"""


def log(message):
    print(message, flush=True)
