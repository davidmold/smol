#!/usr/bin/env python3
"""Trim a GZIPLOGS run's logs back to a checkpoint step, before resuming from it.

    trim_logs.py RUN_DIR STEP [--dry]

With GZIPLOGS each log is one gzip member per checkpoint interval, so the data written
after checkpoint STEP sits in the trailing member(s). For popdy001.dat.gz, species.tsv.gz,
copies.tsv.gz and lives.tsv.gz this keeps the leading members whose rows are all at or
before STEP and drops the rest (including a final member cut off by a crash), so a run
resumed from checkpoint_STEP does not log the same steps twice. The file is truncated in
place; the dropped tail is saved as *.trimmed_tail. A member that straddles STEP is reported and kept whole.

Step column: species.tsv col 2 (born), others col 1 (popdy is comma-separated). Rows within a\nlog are in step order, so only each member's first and last rows are read.
"""
import os, re, sys, zlib

def members(path):
    """Stream the file: [(start, end, first_line, last_line, complete)] per gzip member."""
    out = []
    f = open(path, "rb")
    pos = 0          # file offset of the start of `buf`
    start = 0        # file offset of the current member
    d = zlib.decompressobj(16 + zlib.MAX_WBITS)
    first, tail = None, b""
    buf = f.read(1 << 24)
    while buf:
        try:
            data = d.decompress(buf)
        except zlib.error:
            break
        if data:
            if first is None:
                first = data.split(b"\n", 1)[0]
            tail = (tail + data)[-4096:]
        if d.eof:
            end = pos + len(buf) - len(d.unused_data)
            out.append((start, end, first, tail.rstrip(b"\n").rsplit(b"\n", 1)[-1], True))
            start, pos, buf = end, end, d.unused_data
            d = zlib.decompressobj(16 + zlib.MAX_WBITS)
            first, tail = None, b""
            if not buf:
                pos = end
                buf = f.read(1 << 24)
            continue
        pos += len(buf)
        buf = f.read(1 << 24)
    if first is not None or tail:
        out.append((start, os.path.getsize(path), first or b"", tail.rstrip(b"\n").rsplit(b"\n", 1)[-1], False))
    return out

def step_of(line, col):
    p = re.split(r"[\t, ]+", line.decode(errors="replace").strip())
    try:
        return int(float(p[col]))
    except (IndexError, ValueError):
        return None

def main():
    run, step = sys.argv[1], int(sys.argv[2])
    dry = "--dry" in sys.argv
    for f, col in [("popdy001.dat.gz", 0), ("species.tsv.gz", 1), ("copies.tsv.gz", 0), ("lives.tsv.gz", 0)]:
        path = os.path.join(run, f)
        if not os.path.exists(path):
            continue
        ms = members(path)
        size = os.path.getsize(path)
        keep = 0
        for i, (a0, b0, first, last, ok) in enumerate(ms):
            lo, hi = step_of(first, col), step_of(last, col)
            if ok and hi is not None and hi <= step:
                keep = b0
                continue
            if ok and lo is not None and lo <= step:
                print(f"{f}: member {i} straddles {step} (rows {lo}..{hi}); kept whole")
                keep = b0
                continue
            break
        print(f"{f}: {len(ms)} members, {size} bytes -> keep {keep} bytes")
        if not dry and keep < size:
            with open(path, "rb") as i_, open(path + ".trimmed_tail", "wb") as o:
                i_.seek(keep)
                o.write(i_.read())
            os.truncate(path, keep)

main()
