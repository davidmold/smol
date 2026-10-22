#!/usr/bin/env python3
"""Partnership batch 2 (run on aida): invasion from rarity, bottlenecks, remove-and-rescue,
and accessibility vs introduction in evolving communities. Writes ~/smol/runs/b2/.

    python3 batch2.py [--jobs 9]

All interventions are checkpoint edits (tools/ckpt_edit.py) followed by `smol resume`.
Mutation is off except in the accessibility experiment, which starts from the evolving
strand0 replays (MUTATE 0.0002). Each finished run dir gets log.txt; reruns skip done dirs.
"""
import argparse, os, shutil, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

ap = argparse.ArgumentParser(); ap.add_argument("--jobs", type=int, default=9)
a = ap.parse_args()
H = os.path.expanduser("~/smol"); R = f"{H}/runs"; OUT = f"{R}/b2"
BIN = f"{H}/target/release/smol"; ED = [sys.executable, f"{H}/tools/ckpt_edit.py"]
MK = [sys.executable, f"{H}/experiments/partners/mkconf.py"]
COP = "RERERERE^BQQ$=?>%}"
PAIRS = {"s1": ("S^BO$=J?>G}", "F^O>"), "s2": ("RER$Y$^BQQ$=?>H}", "ERE^AEBP>"),
         "s3": ("FSF^BQ$=?>G}", "SES^AQ"), "g1": ("FUS^YF^BPQ$=?L>}", "SIF^ABQ")}
SEEDS = [18018, 19019, 20020, 21021, 22022, 23023, 24024, 25025, 26026, 27027]
OWN = {18018: ("REER^BQQ$=?>G}", "ERTER^AQQ$T?>G}", [30000, 60000, 90000]),
       20020: ("RECRERDRD^BQQ$=?>%}", "ERRE^AQQ", [20000, 40000, 60000]),
       23023: ("EERRE^BPQ$=?>%}", "RRE^AQQ", [30000, 60000, 90000])}

def sh(cmd, cwd=None):
    subprocess.run(cmd, cwd=cwd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)

def done(d):
    return os.path.exists(f"{d}/log.txt") and os.path.getsize(f"{d}/log.txt") > 0

def resume(d, ck, reseed=None):
    if done(d):
        return
    cmd = [BIN, "resume", ck] + (["--reseed", str(reseed)] if reseed is not None else [])
    with open(f"{d}/log.txt", "w") as f:
        subprocess.run(cmd, cwd=d, stdout=f, stderr=subprocess.STDOUT)

def edit(src_dir, ck, out, ops, rng=1):
    if not os.path.exists(f"{out}/{ck}"):
        sh(ED + [src_dir, ck, out, "--rng-seed", str(rng)] + ops)

def patch(x, y, r=2):
    return ["--region", f"{x - r},{x + r},{y - r},{y + r}", "--clear"]

def pool(jobs):
    with ThreadPoolExecutor(a.jobs) as ex:
        for f in [ex.submit(j) for j in jobs]:
            f.result()

# ---- stage 0: resident copier worlds (mutation off, 50k steps)
def prep_cop(s):
    d = f"{OUT}/copworld/seed{s}"
    if done(d):
        return
    os.makedirs(d, exist_ok=True)
    sh(MK + [f"{d}/run.conf", "--seed", str(s), "--steps", "50000", f"{COP}@62,50"])
    with open(f"{d}/log.txt", "w") as f:
        subprocess.run([BIN, "run", "run.conf"], cwd=d, stdout=f, stderr=subprocess.STDOUT)
pool([lambda s=s: prep_cop(s) for s in SEEDS])

jobs = []
# ---- 1. invasion from rarity (100k steps after introduction)
for p, (A, B) in PAIRS.items():
    for s in SEEDS:
        src = f"{OUT}/copworld/seed{s}"
        for n in (1, 10):  # 1+1 molecules in a 3x3 cleared patch, 10+10 in a 5x5 patch
            d = f"{OUT}/invade/{p}_n{2*n}/seed{s}"; r = 1 if n == 1 else 2
            def j(d=d, src=src, A=A, B=B, n=n, r=r, s=s):
                edit(src, "checkpoint_50000.txt.gz", d, ["--set", "NSTEPS", "150000"] + patch(62, 50, r)
                     + ["--add", A, str(n), "--add", B, str(n)], rng=s)
                resume(d, "checkpoint_50000.txt.gz")
            jobs.append(j)
        # reverse: 20 copiers into a full pair world
        d = f"{OUT}/reverse/{p}/seed{s}"
        def j(d=d, p=p, s=s):
            edit(f"{R}/partners/{p}_AB_m0/seed{s}", "checkpoint_200000.txt.gz", d,
                 ["--set", "NSTEPS", "300000"] + patch(62, 50) + ["--add", COP, "20"], rng=s)
            resume(d, "checkpoint_200000.txt.gz")
        jobs.append(j)

# ---- 2. bottlenecks (cull to 1% and 0.1%, then 100k steps)
for f in (0.01, 0.001):
    for p in PAIRS:
        for s in SEEDS:
            d = f"{OUT}/bottleneck/{p}_f{f}/seed{s}"
            def j(d=d, p=p, s=s, f=f):
                edit(f"{R}/partners/{p}_AB_m0/seed{s}", "checkpoint_200000.txt.gz", d,
                     ["--set", "NSTEPS", "300000", "--cull", str(f)], rng=s)
                resume(d, "checkpoint_200000.txt.gz")
            jobs.append(j)
    for s in SEEDS:
        d = f"{OUT}/bottleneck/cop_f{f}/seed{s}"
        def j(d=d, s=s, f=f):
            edit(f"{OUT}/copworld/seed{s}", "checkpoint_50000.txt.gz", d,
                 ["--set", "NSTEPS", "150000", "--cull", str(f)], rng=s)
            resume(d, "checkpoint_50000.txt.gz")
        jobs.append(j)

# ---- 3. remove one partner (B), add 20 back after a delay, then 100k steps
for p, (A, B) in PAIRS.items():
    for delay in (0, 2000, 5000, 10000, 20000):
        for s in SEEDS[:5]:
            d = f"{OUT}/rescue/{p}_d{delay}/seed{s}"
            def j(d=d, p=p, s=s, B=B, delay=delay):
                src = f"{R}/partners/{p}_AB_m0/seed{s}"
                if delay == 0:
                    edit(src, "checkpoint_200000.txt.gz", d, ["--set", "NSTEPS", "300000", "--remove", B, "--add", B, "20"], rng=s)
                    resume(d, "checkpoint_200000.txt.gz"); return
                mid = d + "/stage1"
                os.makedirs(d, exist_ok=True)
                edit(src, "checkpoint_200000.txt.gz", mid,
                     ["--set", "NSTEPS", str(200000 + delay), "--set", "CKPTEVERY", str(delay), "--remove", B], rng=s)
                resume(mid, "checkpoint_200000.txt.gz")
                ck = f"checkpoint_{200000 + delay}.txt.gz"
                if not os.path.exists(f"{mid}/{ck}"):  # extinct before the delay ran out
                    shutil.copy(f"{mid}/log.txt", f"{d}/log.txt"); return
                edit(mid, ck, d, ["--set", "NSTEPS", "300000", "--set", "CKPTEVERY", "50000", "--add", B, "20"], rng=s)
                resume(d, ck)
            jobs.append(j)

# ---- 4. accessibility vs introduction in evolving communities (MUTATE 0.0002), +150k steps
for s, (A, B, times) in OWN.items():
    src = f"{R}/replay/strand0_seed{s}"
    for t in times:
        ck = f"checkpoint_{t}.txt.gz"
        for k in range(10):  # replays with new random streams: does a partnership form?
            d = f"{OUT}/access/seed{s}_t{t}/replay{k}"
            def j(d=d, src=src, ck=ck, t=t, k=k):
                edit(src, ck, d, ["--set", "NSTEPS", str(t + 150000), "--set", "CKPTKEEPEVERY", "0"])
                resume(d, ck, reseed=1000 + k)
            jobs.append(j)
        for name, (pa, pb) in (("own", (A, B)), ("s1", PAIRS["s1"])):
            for k in range(3):  # the pair introduced at 10+10 molecules in a cleared 5x5 patch
                d = f"{OUT}/access/seed{s}_t{t}/intro_{name}{k}"
                def j(d=d, src=src, ck=ck, t=t, k=k, pa=pa, pb=pb):
                    edit(src, ck, d, ["--set", "NSTEPS", str(t + 150000), "--set", "CKPTKEEPEVERY", "0"] + patch(30 + 30 * k, 50)
                         + ["--add", pa, "10", "--add", pb, "10"], rng=k + 1)
                    resume(d, ck, reseed=2000 + k)
                jobs.append(j)

print(f"{len(jobs)} jobs", flush=True)
pool(jobs)
open(f"{OUT}/DONE", "w").write("done\n")
