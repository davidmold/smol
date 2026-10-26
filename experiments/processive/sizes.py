#!/usr/bin/env python3
"""Introduction size: how many processive molecules does it take to establish? (run on aida)

    sizes.py SEQS.json [--jobs 8] [--steps 200000]

SEQS.json: {"name": seq, ...}: big-world processive genotypes and the processive variants seen
briefly in small worlds (findings 14).
Recipients: small classic communities from pr0x at 2.5M (COMM below). Every treatment clears the
same 5x5 patch at (60-64, 48-52); then N molecules of the sequence are added inside a box of
that size centred on (62, 50):
    N=1   (62,50)
    N=2   (62,50),(63,50)            adjacent pair
    N=5   5 of the 9 cells 61-63 x 49-51
    N=20  20 of the 25 cells of the patch (the earlier transplant protocol)
    NONE  patch cleared, nothing added (the matched clearing control)
Replicates per recipient x sequence: N=1 and N=2 x12, N=5 x6, N=20 x3; NONE x3. Mutation on,
OFFSPRINGLOG + GZIPLOGS + SPECIESFORGET, and only the live species are copied in
(ckpt_edit.py --live-species). Output: ~/smol/runs/sizes/<recipient>/<seq>/n<N>/rep<k>
"""
import argparse, json, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

ap = argparse.ArgumentParser(); ap.add_argument("seqs")
ap.add_argument("--jobs", type=int, default=8); ap.add_argument("--steps", type=int, default=200000)
a = ap.parse_args()
S = json.load(open(a.seqs))
H = os.path.expanduser("~/smol"); OUT = f"{H}/runs/sizes"
BIN = f"{H}/target/release/smol"; ED = [sys.executable, f"{H}/tools/ckpt_edit.py"]
COMM = ["seed36036", "seed21021"]
BOX = {1: "62,62,50,50", 2: "62,63,50,50", 5: "61,63,49,51", 20: "60,64,48,52"}
REPS = {1: 12, 2: 12, 5: 6, 20: 3}

def done(d):
    return os.path.exists(f"{d}/log.txt") and os.path.getsize(f"{d}/log.txt") > 0

def intro(comm, name, seq, n, k):
    d = f"{OUT}/{comm}/{name}/n{n}/rep{k}"
    if done(d):
        return
    src = f"{H}/runs/pr0x/{comm}"; ck = "checkpoint_2500000.txt.gz"
    ops = ["--live-species", "--set", "NSTEPS", str(2500000 + a.steps), "--set", "GZIPLOGS", "1",
           "--set", "SPECIESFORGET", "1", "--set", "OFFSPRINGLOG", "1", "--set", "CKPTEVERY", "50000",
           "--set", "CKPTKEEPEVERY", "100000", "--region", "60,64,48,52", "--clear"]
    if seq:
        ops += ["--region", BOX[n], "--add", seq, str(n)]
    subprocess.run(ED + [src, ck, d, "--rng-seed", str(k + 1)] + ops, check=True, stdout=subprocess.DEVNULL)
    with open(f"{d}/log.txt", "w") as f:
        subprocess.run([BIN, "resume", ck, "--reseed", str(9300 + k)], cwd=d, stdout=f, stderr=subprocess.STDOUT)

jobs = []
for c in COMM:
    for k in range(3):
        jobs.append(lambda c=c, k=k: intro(c, "NONE", None, 0, k))
    for n in (20, 5, 2, 1):
        for k in range(REPS[n]):
            for name, seq in S.items():
                jobs.append(lambda c=c, nm=name, s=seq, n=n, k=k: intro(c, nm, s, n, k))
with ThreadPoolExecutor(a.jobs) as ex:
    for f in [ex.submit(j) for j in jobs]:
        f.result()
open(f"{OUT}/DONE", "w").write("done\n")
