#!/usr/bin/env python3
"""Edit a smol checkpoint to set up an intervention, then resume it in a new run directory.

    ckpt_edit.py RUN_DIR CKPT OUT_DIR [ops...] [--rng-seed N]

Ops, applied in the order given:
  --remove SEQ           remove every molecule whose sequence is exactly SEQ
  --remove-species ID    remove every molecule of that species id
  --cull F               keep each molecule with probability F (a bottleneck)
  --add SEQ N            put N new free molecules of SEQ in random empty cells
  --region X0,X1,Y0,Y1   restrict the following --add ops to that box (inclusive; storage x)
  --set KEY VALUE        set a config key in the checkpoint (e.g. MUTATE 0, NSTEPS 400000)

A molecule whose partner is removed is freed (its reaction ends, any half-built copy is
dropped). OUT_DIR gets the edited checkpoint and a copy of RUN_DIR/species.tsv (plus any new
species); then:  cd OUT_DIR && smol resume <checkpoint> [--reseed S]
Works on plain or gzipped species.tsv.
"""
import argparse, gzip, os, random, sys

ap = argparse.ArgumentParser()
ap.add_argument("run"); ap.add_argument("ckpt"); ap.add_argument("out")
ap.add_argument("--rng-seed", type=int, default=1)
a, ops = ap.parse_known_args()
rng = random.Random(a.rng_seed)

def ropen(p):
    if os.path.exists(p):
        return open(p)
    return gzip.open(p + ".gz", "rt")

lines = gzip.open(os.path.join(a.run, a.ckpt), "rt").read().split("\n")
assert lines[0] == "SMOLCKPT 1"
step = int(lines[1].split()[1]); rngl = lines[2]
nextsp = int(lines[3].split()[1]); nconf = int(lines[4].split()[1])
conf = lines[5:5 + nconf]
i = 5 + nconf
slots = int(lines[i].split()[1]); free = [int(v) for v in lines[i + 1].split()[2:]]
nm = int(lines[i + 2].split()[1])
mols = {}  # idx -> [x, y, species, seq, bondtokens]
for l in lines[i + 3:i + 3 + nm]:
    p = l.split()
    mols[int(p[0])] = [int(p[1]), int(p[2]), int(p[3]), p[4], p[5:]]

def cval(k, d):
    for l in conf:
        t = l.split()
        if t and t[0] == k:
            return t[1]
    return d
gx, gy, depth = int(float(cval("GRIDX", 125))), int(float(cval("GRIDY", 100))), int(float(cval("DEPTH", 1)))
W = gx * depth

species = {}
spfile = []
for l in ropen(os.path.join(a.run, "species.tsv")):
    spfile.append(l)
    p = l.rstrip("\n").split("\t")
    if p[0].isdigit():
        species.setdefault(p[-1], int(p[0]))

def partner(idx):
    b = mols[idx][4]
    return int(b[1]) if b[0] in ("A", "P") else None

def remove(idx):
    q = partner(idx)
    if q is not None and q in mols:
        mols[q][4] = ["F"]
    del mols[idx]
    free.append(idx)

region = None; added = 0; removed = 0
k = 0
while k < len(ops):
    op = ops[k]
    if op == "--remove":
        for idx in [j for j, m in mols.items() if m[3].replace(".", "") == ops[k + 1]]:
            if idx in mols: remove(idx); removed += 1
        k += 2
    elif op == "--remove-species":
        sp = int(ops[k + 1])
        for idx in [j for j, m in mols.items() if m[2] == sp]:
            if idx in mols: remove(idx); removed += 1
        k += 2
    elif op == "--cull":
        f = float(ops[k + 1])
        for idx in list(mols):
            if idx in mols and rng.random() >= f:
                remove(idx); removed += 1
        k += 2
    elif op == "--region":
        region = [int(v) for v in ops[k + 1].split(",")]; k += 2
    elif op == "--add":
        seq, n = ops[k + 1], int(ops[k + 2])
        if seq not in species:
            species[seq] = nextsp
            spfile.append(f"{nextsp}\t{step}\t-1\t-1\t{seq}\n"); nextsp += 1
        occ = {(m[0], m[1]) for m in mols.values()}
        x0, x1, y0, y1 = region if region else (0, W - 1, 0, gy - 1)
        empty = [(x, y) for x in range(x0, x1 + 1) for y in range(y0, y1 + 1) if (x, y) not in occ]
        if len(empty) < n:
            sys.exit(f"only {len(empty)} empty cells for --add {seq} {n}")
        for (x, y) in rng.sample(empty, n):
            idx = free.pop(0) if free else slots
            if idx == slots:
                slots += 1
            mols[idx] = [x, y, species[seq], seq, ["F"]]
            added += 1
        k += 3
    elif op == "--set":
        key, val = ops[k + 1], ops[k + 2]
        hit = False
        for j, l in enumerate(conf):
            t = l.split()
            if t and t[0] == key:
                conf[j] = f"{key} {val}"; hit = True
        if not hit:
            conf.insert(0, f"{key} {val}")
        k += 3
    else:
        sys.exit(f"unknown op {op}")

os.makedirs(a.out, exist_ok=True)
name = os.path.basename(a.ckpt)
with gzip.open(os.path.join(a.out, name), "wt") as f:
    f.write(f"SMOLCKPT 1\nSTEP {step}\n{rngl}\nNEXTSPECIES {nextsp}\nCONFIG {len(conf)}\n")
    for l in conf:
        f.write(l + "\n")
    f.write(f"SLOTS {slots}\nFREE {len(free)}" + "".join(f" {v}" for v in free) + "\n")
    f.write(f"MOLS {len(mols)}\n")
    for idx in sorted(mols):
        x, y, sp, seq, b = mols[idx]
        f.write(f"{idx} {x} {y} {sp} {seq} {' '.join(b)}\n")
with open(os.path.join(a.out, "species.tsv"), "w") as f:
    f.writelines(spfile)
print(f"step {step}: removed {removed}, added {added}, molecules now {len(mols)} -> {a.out}/{name}")
