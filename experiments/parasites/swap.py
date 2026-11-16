#!/usr/bin/env python3
"""swap.py RECIP_SRC RECIP_STEP DONOR_SRC DONOR_STEP OUT SEED NSTEPS [CLASS]
Replace every molecule of CLASS (default R, replicators; P for parasites) in RECIP_SRC/classes.tsv in the recipient world with
a free molecule drawn, by abundance, from the donor world's replicators (class R in
DONOR_SRC/classes.tsv). Parasites and other molecules stay where they are; the partners of
replaced molecules are freed. With DONOR = RECIP this is the control: same disturbance, same
kind of replicators."""
import sys, gzip, random, subprocess, os, collections
rs, rstep, ds, dstep, out, seed, n = sys.argv[1], int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), sys.argv[5], int(sys.argv[6]), int(sys.argv[7])
K = sys.argv[8] if len(sys.argv) > 8 else "R"
rng = random.Random(seed)
def classes(d):
    return {int(l.split()[0]): l.split()[1] for l in open(f"{d}/classes.tsv")}
def seqs(d):
    s = {}
    for l in open(f"{d}/species.tsv"):
        p = l.rstrip("\n").split("\t")
        if p[0].isdigit(): s.setdefault(int(p[0]), p[-1])
    return s
def mols(d, step):
    with gzip.open(f"{d}/checkpoint_{step}.txt.gz", "rt") as f:
        for l in f:
            if l.startswith("MOLS"): break
        for l in f:
            p = l.split(); yield int(p[1]), int(p[2]), int(p[3])
rc, dc, dseq = classes(rs), classes(ds), seqs(ds)
pool = [dseq[sp] for x, y, sp in mols(ds, dstep) if dc.get(sp) == K]
cells = [(x, y) for x, y, sp in mols(rs, rstep) if rc.get(sp) == K]
ops = []
for (x, y) in cells:
    ops += ["--put", f"{x},{y}", rng.choice(pool)]
os.makedirs(out, exist_ok=True)
cmd = ["python3", os.path.expanduser("~/smol-src/tools/ckpt_edit.py"), rs, f"checkpoint_{rstep}.txt.gz", out] + ops + \
      ["--set", "NSTEPS", str(rstep + n), "--set", "POPDYEVERY", "1000", "--rng-seed", str(seed)]
r = subprocess.run(cmd, capture_output=True, text=True)
open(f"{out}/edit.txt", "w").write(f"replaced {len(cells)} {K} cells from a pool of {len(pool)}\n" + r.stdout + r.stderr)
sys.exit(r.returncode)
