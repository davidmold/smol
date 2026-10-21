#!/usr/bin/env python3
"""Configs for the partnership isolation experiments.

    mkconf.py OUT.conf --seed S --mutate M --steps N  BLOCK [BLOCK ...]

BLOCK = SEQ_A[+SEQ_B]@X,Y[,SIZE]  : an SIZE x SIZE block (default 11) centred at X,Y; with
two sequences they alternate in a checkerboard. Copy-strand chemistry, 125x100, OFFSPRINGLOG.
"""
import argparse
ap = argparse.ArgumentParser()
ap.add_argument("out"); ap.add_argument("blocks", nargs="+")
ap.add_argument("--seed", type=int, default=18018); ap.add_argument("--mutate", default="0")
ap.add_argument("--steps", type=int, default=200000); ap.add_argument("--grid", default="125x100")
ap.add_argument("--extra", action="append", default=[])
a = ap.parse_args()
gx, gy = map(int, a.grid.split("x"))
L = [f"GRIDX {gx}", f"GRIDY {gy}", "DECAY 0.0005", f"MUTATE {a.mutate}", f"RANDSEED {a.seed}",
     f"NSTEPS {a.steps}", "STRANDCOPY 1", "OFFSPRINGLOG 1", "POPDYEVERY 1000", "SNAPEVERY 2000",
     "COPYEVERY 1000", "CKPTEVERY 50000", "KEEPCKPT 2", "CKPTKEEPEVERY 50000"] + a.extra
taken = set()
for b in a.blocks:
    seqs, where = b.rsplit("@", 1)
    parts = [int(v) for v in where.split(",")]
    cx, cy = parts[0], parts[1]; size = parts[2] if len(parts) > 2 else 11
    seqs = seqs.split("+")
    cells = {s: [] for s in seqs}
    for i in range(size):
        for j in range(size):
            x, y = (cx - size // 2 + i) % gx, (cy - size // 2 + j) % gy
            if (x, y) in taken:
                continue
            taken.add((x, y))
            cells[seqs[(i + j) % len(seqs)]].append((x, y))
    for s, cs in cells.items():
        L.append(f"AGENT  {s} {len(cs)} Q")
        L += [f"GRIDPOS {x} {y}" for x, y in cs]
open(a.out, "w").write("\n".join(L) + "\n")
