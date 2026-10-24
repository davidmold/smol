#!/usr/bin/env python3
"""Is the evolved community still mutually dependent? Removal tests on the descendants.

    dep.py [--jobs 10]       (on aida; writes ~/smol/runs/dep/)

For every mutation-on pair-alone run of batch 1 (s1, s2, s3, g1 x 10 seeds) that is alive at
200k, from its 200k checkpoint (mutation stays on), 50k more steps under four treatments:
  ctrl     unchanged
  noB      every species of the B lineage removed
  noA      every species of the A lineage removed
  rescueB  B lineage removed, then 20 molecules of the most abundant B-lineage species that
           made placed copies (a functional evolved partner) added back at random
Lineages come from template ancestry (as in lineage.py).
"""
import argparse, gzip, json, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

ap = argparse.ArgumentParser(); ap.add_argument("--jobs", type=int, default=10)
a = ap.parse_args()
H = os.path.expanduser("~/smol"); P = f"{H}/runs/partners"; OUT = f"{H}/runs/dep"
BIN = f"{H}/target/release/smol"; ED = [sys.executable, f"{H}/tools/ckpt_edit.py"]
PAIRS = {"s1": ("S^BO$=J?>G}", "F^O>"), "s2": ("RER$Y$^BQQ$=?>H}", "ERE^AEBP>"),
         "s3": ("FSF^BQ$=?>G}", "SES^AQ"), "g1": ("FUS^YF^BPQ$=?L>}", "SIF^ABQ")}

def ropen(p):
    return open(p) if os.path.exists(p) else gzip.open(p + ".gz", "rt")

def lineages(run, founders):
    seq, tpl = {}, {}
    for l in ropen(f"{run}/species.tsv"):
        p = l.rstrip("\n").split("\t")
        if p[0].isdigit():
            seq[int(p[0])] = p[-1]; tpl[int(p[0])] = int(p[3])
    roots = {v: k for k, v in founders.items()}; memo = {}
    def lin(s):
        path = []
        while s not in memo:
            if seq.get(s) in roots and tpl.get(s, -1) == -1:
                memo[s] = roots[seq[s]]; break
            path.append(s); p = tpl.get(s, -1)
            if p == -1 or p == s:
                memo[s] = "unknown"; break
            s = p
        for x in path:
            memo[x] = memo[s]
        return memo[s]
    return seq, lin

def ckpt_species(run, ck):
    ids = set(); mols = False
    for l in gzip.open(f"{run}/{ck}", "rt"):
        if l.startswith("MOLS"):
            mols = True; continue
        if mols:
            p = l.split()
            if len(p) > 4:
                ids.add(int(p[3]))
    return ids

def job(p, seed):
    run = f"{P}/{p}_AB_m0.0002/seed{seed}"; ck = "checkpoint_200000.txt.gz"
    if "EXTINCT" in open(f"{run}/log.txt").read() or not os.path.exists(f"{run}/{ck}"):
        return
    A, B = PAIRS[p]
    seq, lin = lineages(run, {"A": A, "B": B})
    ids = ckpt_species(run, ck)
    Bids = [s for s in ids if lin(s) == "B"]; Aids = [s for s in ids if lin(s) == "A"]
    roles = json.load(open(f"{P}/lin/{p}_AB_seed{seed}.json"))["roles"]
    fB = [r["seq"] for r in roles if r["lin"] == "B" and r["made"] > 0]
    treat = {"ctrl": [], "noB": sum((["--remove-species", str(s)] for s in Bids), []),
             "noA": sum((["--remove-species", str(s)] for s in Aids), [])}
    if fB:
        treat["rescueB"] = treat["noB"] + ["--add", fB[0], "20"]
    for t, ops in treat.items():
        d = f"{OUT}/{p}/{t}/seed{seed}"
        if os.path.exists(f"{d}/log.txt"):
            continue
        subprocess.run(ED + [run, ck, d, "--set", "NSTEPS", "250000", "--set", "CKPTKEEPEVERY", "0", "--rng-seed", str(seed)] + ops,
                       check=True, stdout=subprocess.DEVNULL)
        json.dump({"A_species": len(Aids), "B_species": len(Bids), "rescue_seq": fB[0] if fB else None},
                  open(f"{d}/treat.json", "w"))
        with open(f"{d}/log.txt", "w") as f:
            subprocess.run([BIN, "resume", ck], cwd=d, stdout=f, stderr=subprocess.STDOUT)

with ThreadPoolExecutor(a.jobs) as ex:
    fs = [ex.submit(job, p, s) for p in PAIRS for s in range(18018, 28028, 1001)]
    for f in fs:
        f.result()
open(f"{OUT}/DONE", "w").write("done\n")
