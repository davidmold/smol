#!/usr/bin/env python3
"""Replays of pr0xx/seed21021 near the natural emergence of its processive lineage (run on aida).

    replays.py [--jobs 7] [--steps 200000] [--reps 24]

From the 9.9M checkpoint (before S21's lineage arose at 9.96-9.98M) and the 10.0M checkpoint
(lineage present, lost by 10.02M in the original run), resume with a different random stream
(--reseed) and nothing else changed, for 200k steps. How often does a processive population
establish when the same community is replayed? This is the no-introduction control for
natural.py, two of whose E1-introduction runs from 9.9M ended with an established native
processive population (S21-like, not E1-derived).
Mutation on, OFFSPRINGLOG, GZIPLOGS, SPECIESFORGET; species from natural.py's live-species
copies (runs/natural/_src). Output: ~/smol/runs/replays/<s21_9.9M|s21_10.0M>/rep<k>
"""
import argparse, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

ap = argparse.ArgumentParser()
ap.add_argument("--jobs", type=int, default=7); ap.add_argument("--steps", type=int, default=200000)
ap.add_argument("--reps", type=int, default=24)
a = ap.parse_args()
H = os.path.expanduser("~/smol"); OUT = f"{H}/runs/replays"
BIN = f"{H}/target/release/smol"; ED = [sys.executable, f"{H}/tools/ckpt_edit.py"]
REC = {"s21_9.9M": 9900000, "s21_10.0M": 10000000}

def replay(rec, k):
    d = f"{OUT}/{rec}/rep{k}"
    if os.path.exists(f"{d}/log.txt") and os.path.getsize(f"{d}/log.txt") > 0:
        return
    step = REC[rec]; ck = f"checkpoint_{step}.txt.gz"
    ops = ["--live-species", "--set", "NSTEPS", str(step + a.steps), "--set", "GZIPLOGS", "1",
           "--set", "SPECIESFORGET", "1", "--set", "OFFSPRINGLOG", "1", "--set", "CKPTEVERY", "50000",
           "--set", "KEEPCKPT", "1", "--set", "CKPTKEEPEVERY", "0"]
    os.makedirs(d, exist_ok=True)
    subprocess.run(ED + [f"{H}/runs/natural/_src/{rec}", ck, d] + ops, check=True, stdout=subprocess.DEVNULL)
    with open(f"{d}/log.txt", "w") as f:
        subprocess.run([BIN, "resume", ck, "--reseed", str(9700 + k)], cwd=d, stdout=f, stderr=subprocess.STDOUT)

jobs = [lambda r=r, k=k: replay(r, k) for k in range(a.reps) for r in REC]
with ThreadPoolExecutor(a.jobs) as ex:
    for f in [ex.submit(j) for j in jobs]:
        f.result()
open(f"{OUT}/DONE", "w").write("done\n")
