#!/usr/bin/env python3
"""Does spatially clustered disturbance help spontaneous innovations establish? (run on aida)

    disturb.py [--jobs 10] [--reps 40] [--steps 200000] [--until 150000]

Start: pr0xx/seed21021 at 9.9M (the live-species copy in runs/natural/_src/s21_9.9M), where a
spontaneous processive transition happened in 1 of 24 clean replays. Nothing is introduced and
the chemistry is unchanged. Treatments, each with --reps reseeded runs:
    base       existing dynamics
    clustered  every 1,000 steps remove every molecule in one random 7x7 patch (DISTURB_MODE 1)
    dispersed  every 1,000 steps remove whatever is in 49 random cells (DISTURB_MODE 2): the same
               expected mortality and freed space, without spatial clustering
Disturbance is blind to genotype and to any detected novelty, runs for the first --until steps
and then stops; every run goes on to --steps. Rep k uses the same random stream (--reseed
9800+k) and the same disturbance seed (1000+k) in every treatment. Checkpoints every 50k are
kept (for continuing populations that establish). Output: ~/smol/runs/disturb/<treatment>/rep<k>
"""
import argparse, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

ap = argparse.ArgumentParser()
ap.add_argument("--jobs", type=int, default=10); ap.add_argument("--reps", type=int, default=40)
ap.add_argument("--steps", type=int, default=200000); ap.add_argument("--until", type=int, default=150000)
a = ap.parse_args()
H = os.path.expanduser("~/smol"); OUT = f"{H}/runs/disturb"; SRC = f"{H}/runs/natural/_src/s21_9.9M"
BIN = f"{H}/target/release/smol"; ED = [sys.executable, f"{H}/tools/ckpt_edit.py"]
STEP = 9900000; CK = f"checkpoint_{STEP}.txt.gz"
TREAT = {"base": None, "clustered": 1, "dispersed": 2}

def run(treat, k):
    d = f"{OUT}/{treat}/rep{k}"
    if os.path.exists(f"{d}/log.txt") and os.path.getsize(f"{d}/log.txt") > 0:
        return
    ops = ["--live-species", "--set", "NSTEPS", str(STEP + a.steps), "--set", "GZIPLOGS", "1",
           "--set", "SPECIESFORGET", "1", "--set", "OFFSPRINGLOG", "1", "--set", "CKPTEVERY", "50000",
           "--set", "KEEPCKPT", "1", "--set", "CKPTKEEPEVERY", "50000"]
    mode = TREAT[treat]
    if mode:
        ops += ["--set", "DISTURB_EVERY", "1000", "--set", "DISTURB_SIZE", "7", "--set", "DISTURB_MODE", str(mode),
                "--set", "DISTURB_UNTIL", str(STEP + a.until), "--set", "DISTURB_SEED", str(1000 + k)]
    os.makedirs(d, exist_ok=True)
    subprocess.run(ED + [SRC, CK, d] + ops, check=True, stdout=subprocess.DEVNULL)
    with open(f"{d}/log.txt", "w") as f:
        subprocess.run([BIN, "resume", CK, "--reseed", str(9800 + k)], cwd=d, stdout=f, stderr=subprocess.STDOUT)

jobs = [lambda t=t, k=k: run(t, k) for k in range(a.reps) for t in TREAT]
with ThreadPoolExecutor(a.jobs) as ex:
    for f in [ex.submit(j) for j in jobs]:
        f.result()
open(f"{OUT}/DONE", "w").write("done\n")
