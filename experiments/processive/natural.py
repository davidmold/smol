#!/usr/bin/env python3
"""Naturalistic introductions: one (or two adjacent) processive molecules replacing residents,
no cleared patch, at a different random location in each replicate, with the newcomers'
lineage traced (run on aida).

    natural.py SEQS.json [--jobs 7] [--steps 200000] [--reps 16]

Recipients (checkpoint, community):
    s36_2.5M    pr0x/seed36036 at 2.5M   (permissive in sizes.py)
    s21_2.5M    pr0x/seed21021 at 2.5M   (resistant in sizes.py; ancestor of the next two)
    s21_9.9M    pr0xx/seed21021 at 9.9M  (just before S21's lineage arose there, 9.96-9.98M)
    s21_10.0M   pr0xx/seed21021 at 10.0M (S21's lineage present, about to be lost)
Treatments: N=1 replaces the molecule in one randomly chosen occupied cell (its partner, if it
had one, is freed); N=2 also replaces the cell to its right (x+1, wrapping). Nothing else is
touched. Sequences: SEQS.json (E1 = big-world processive, S21 = seed21021's own variant).
Each run: mutation on, OFFSPRINGLOG, GZIPLOGS, SPECIESFORGET, live species only, and
TRACESEQ = the introduced sequence (trace.tsv: binds, reaction ends, products, lost products,
deaths of the newcomers and their descendants; TRACEMAX 100000 lines).
Output: ~/smol/runs/natural/<recipient>/<seq>/n<N>/rep<k>; intro.txt records the cells used.
Both sequences and both N share a location in a given (recipient, rep), so treatments are paired.
"""
import argparse, gzip, json, os, random, subprocess, sys, zlib
from concurrent.futures import ThreadPoolExecutor

ap = argparse.ArgumentParser(); ap.add_argument("seqs")
ap.add_argument("--jobs", type=int, default=7); ap.add_argument("--steps", type=int, default=200000)
ap.add_argument("--reps", type=int, default=16)
a = ap.parse_args()
S = json.load(open(a.seqs))
H = os.path.expanduser("~/smol"); OUT = f"{H}/runs/natural"
BIN = f"{H}/target/release/smol"; ED = [sys.executable, f"{H}/tools/ckpt_edit.py"]
REC = {"s36_2.5M": ("pr0x/seed36036", 2500000), "s21_2.5M": ("pr0x/seed21021", 2500000),
       "s21_9.9M": ("pr0xx/seed21021", 9900000), "s21_10.0M": ("pr0xx/seed21021", 10000000)}

def prep(rec):
    """A one-off copy of each recipient checkpoint with only its live species, so the runs don't
    each stream the community's whole species history."""
    src, step = REC[rec]; d = f"{OUT}/_src/{rec}"
    if not os.path.exists(f"{d}/checkpoint_{step}.txt.gz"):
        subprocess.run(ED + [f"{H}/runs/{src}", f"checkpoint_{step}.txt.gz", d, "--live-species"], check=True,
                       stdout=subprocess.DEVNULL)
    return d

def occupied(src, step):
    """(gx, occupied cells) of a checkpoint."""
    L = gzip.open(f"{H}/runs/{src}/checkpoint_{step}.txt.gz", "rt").read().split("\n")
    nconf = int(L[4].split()[1]); conf = L[5:5 + nconf]
    gx = next(int(float(l.split()[1])) for l in conf if l.split()[:1] == ["GRIDX"])
    i = 5 + nconf; nm = int(L[i + 2].split()[1])
    return gx, sorted({(int(p.split()[1]), int(p.split()[2])) for p in L[i + 3:i + 3 + nm]})

OCC = {r: occupied(*REC[r]) for r in REC}
with ThreadPoolExecutor(4) as ex:
    SRC = dict(zip(REC, ex.map(prep, REC)))

def done(d):
    return os.path.exists(f"{d}/log.txt") and os.path.getsize(f"{d}/log.txt") > 0

def intro(rec, name, seq, n, k):
    d = f"{OUT}/{rec}/{name}/n{n}/rep{k}"
    if done(d):
        return
    src, step = REC[rec]; ck = f"checkpoint_{step}.txt.gz"
    gx, occ = OCC[rec]
    # one location per (recipient, rep), shared by both sequences and both N, so treatments are paired
    x, y = random.Random(zlib.crc32(f"{rec}/{k}".encode())).choice(occ)
    cells = [(x, y)] + ([((x + 1) % gx, y)] if n == 2 else [])
    ops = ["--live-species", "--set", "NSTEPS", str(step + a.steps), "--set", "GZIPLOGS", "1",
           "--set", "SPECIESFORGET", "1", "--set", "OFFSPRINGLOG", "1", "--set", "CKPTEVERY", "50000",
           "--set", "KEEPCKPT", "1", "--set", "CKPTKEEPEVERY", "0",
           "--set", "TRACESEQ", seq, "--set", "TRACEMAX", "100000",
           "--set", "TRACEAT", ";".join(f"{cx},{cy}" for cx, cy in cells)]
    for cx, cy in cells:
        ops += ["--put", f"{cx},{cy}", seq]
    os.makedirs(d, exist_ok=True)
    open(f"{d}/intro.txt", "w").write(f"{rec} {src} {step} {name} n={n} cells={cells}\n")
    subprocess.run(ED + [SRC[rec], ck, d, "--rng-seed", str(k + 1)] + ops, check=True, stdout=subprocess.DEVNULL)
    with open(f"{d}/log.txt", "w") as f:
        subprocess.run([BIN, "resume", ck, "--reseed", str(9500 + k)], cwd=d, stdout=f, stderr=subprocess.STDOUT)

jobs = [lambda r=r, nm=nm, s=s, n=n, k=k: intro(r, nm, s, n, k)
        for k in range(a.reps) for r in REC for n in (1, 2) for nm, s in S.items()]
with ThreadPoolExecutor(a.jobs) as ex:
    for f in [ex.submit(j) for j in jobs]:
        f.result()
open(f"{OUT}/DONE", "w").write("done\n")
