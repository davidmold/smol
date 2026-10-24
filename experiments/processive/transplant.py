#!/usr/bin/env python3
"""Transplant processive replicators into small classic worlds (run on aida).

    transplant.py SEQS.json [--jobs 9] [--steps 300000]

SEQS.json: {"E1": seq, "E2": seq, "M1": seq, "M2": seq, "ORD": seq}
  E = early processive (from 4.5-5M of big/classic_400x300), M = mature (9.5M),
  ORD = an ordinary (one copy per reaction) replicator from before the transition (control).
Recipients:
  empty   a fresh 125x100 world containing only the 20 introduced molecules
  pr0x/seedS  six small classic communities evolved for 2.5M steps (from their 2.5M checkpoint)
Introduction: clear a 5x5 patch at (62,50) and add 20 molecules of the sequence (adjacent, so
failure is not failure to meet a partner); NONE = patch cleared, nothing added (resident-only
control). Mutation on (0.0002), OFFSPRINGLOG + GZIPLOGS, 3 replicates (different random
streams) per recipient x treatment. Output: ~/smol/runs/transplant/<recipient>/<treatment>/rep<k>
"""
import argparse, json, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

ap = argparse.ArgumentParser(); ap.add_argument("seqs")
ap.add_argument("--jobs", type=int, default=9); ap.add_argument("--steps", type=int, default=300000)
a = ap.parse_args()
S = json.load(open(a.seqs))
H = os.path.expanduser("~/smol"); OUT = f"{H}/runs/transplant"
BIN = f"{H}/target/release/smol"; ED = [sys.executable, f"{H}/tools/ckpt_edit.py"]
COMM = ["seed21021", "seed29029", "seed36036", "seed40040", "seed47047", "seed56056"]

def done(d):
    return os.path.exists(f"{d}/log.txt") and os.path.getsize(f"{d}/log.txt") > 0

def empty(name, seq, k):
    d = f"{OUT}/empty/{name}/rep{k}"
    if done(d):
        return
    os.makedirs(d, exist_ok=True)
    cells = [(60 + i, 48 + j) for i in range(5) for j in range(5)][:20]
    L = ["GRIDX 125", "GRIDY 100", "DECAY 0.0005", "MUTATE 0.0002", f"RANDSEED {7001 + k}",
         f"NSTEPS {a.steps}", "OFFSPRINGLOG 1", "GZIPLOGS 1", "POPDYEVERY 1000", "SNAPEVERY 5000",
         "CKPTEVERY 50000", "KEEPCKPT 2", "CKPTKEEPEVERY 100000", f"AGENT  {seq} 20 Q"]
    L += [f"GRIDPOS {x} {y}" for x, y in cells]
    open(f"{d}/run.conf", "w").write("\n".join(L) + "\n")
    with open(f"{d}/log.txt", "w") as f:
        subprocess.run([BIN, "run", "run.conf"], cwd=d, stdout=f, stderr=subprocess.STDOUT)

def intro(comm, name, seq, k):
    d = f"{OUT}/{comm}/{name}/rep{k}"
    if done(d):
        return
    src = f"{H}/runs/pr0x/{comm}"; ck = "checkpoint_2500000.txt.gz"
    ops = ["--set", "NSTEPS", str(2500000 + a.steps), "--set", "GZIPLOGS", "1", "--set", "CKPTEVERY", "50000",
           "--set", "CKPTKEEPEVERY", "100000", "--region", "60,64,48,52", "--clear"]
    if seq:
        ops += ["--add", seq, "20"]
    subprocess.run(ED + [src, ck, d, "--rng-seed", str(k + 1)] + ops, check=True, stdout=subprocess.DEVNULL)
    with open(f"{d}/log.txt", "w") as f:
        subprocess.run([BIN, "resume", ck, "--reseed", str(9100 + k)], cwd=d, stdout=f, stderr=subprocess.STDOUT)

jobs = []
for k in range(3):
    for name, seq in S.items():
        jobs.append(lambda n=name, s=seq, k=k: empty(n, s, k))
    for c in COMM:
        for name, seq in list(S.items()) + [("NONE", None)]:
            jobs.append(lambda c=c, n=name, s=seq, k=k: intro(c, n, s, k))
with ThreadPoolExecutor(a.jobs) as ex:
    for f in [ex.submit(j) for j in jobs]:
        f.result()
open(f"{OUT}/DONE", "w").write("done\n")
