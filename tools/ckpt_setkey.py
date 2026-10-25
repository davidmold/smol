#!/usr/bin/env python3
"""Set (or add) a config key in a checkpoint, in place of ckpt_edit.py when nothing else changes.

    ckpt_setkey.py IN.txt.gz OUT.txt.gz KEY VALUE

Resuming takes its config from the checkpoint, so this is how a switch such as
SPECIESFORGET 1 is turned on for the rest of a run. Writes OUT next to IN; resume from OUT
in the same run directory so the logs carry on.
"""
import gzip, sys

src, dst, key, val = sys.argv[1:5]
with gzip.open(src, "rt") as f:
    lines = f.read().split("\n")
ci = next(i for i, l in enumerate(lines) if l.startswith("CONFIG "))
n = int(lines[ci][7:])
conf = lines[ci + 1:ci + 1 + n]
hit = False
for j, l in enumerate(conf):
    p = l.split()
    if p and p[0] == key:
        conf[j] = f"{key} {val}"
        hit = True
if not hit:
    conf.append(f"{key} {val}")
lines[ci:ci + 1 + n] = [f"CONFIG {len(conf)}"] + conf
with gzip.open(dst, "wt", compresslevel=3) as f:
    f.write("\n".join(lines))
print(("set" if hit else "added"), key, val, "->", dst)
