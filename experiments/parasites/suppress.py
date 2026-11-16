#!/usr/bin/env python3
"""suppress.py OUT MODE SEED NSTEPS [K]

A whole run from the seed replicators (the paper's "box" config, old rules) in K-step chunks.
After each chunk every species present is classed from the run's copy log so far:
  copier  ever made an exact copy (of itself or a partner)
  P       was copied exactly by someone, never made a copy itself, and has existed for >= K steps
MODE:
  none   no edit (the same restart every chunk, as a control for the restarts)
  P      remove every molecule of a P species   (parasites suppressed throughout)
  cull   remove the same number of molecules at random (control for the removal itself)
Every chunk logs: step, molecules, P molecules, molecules removed  -> suppress.tsv
Checkpoints every 100k are kept (CKPTKEEPEVERY) for later analysis."""
import sys, os, gzip, subprocess, collections, glob

out, mode, seed, nsteps = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
K = int(sys.argv[5]) if len(sys.argv) > 5 else 5000
H = os.path.expanduser("~")
SMOL = f"{H}/smol-src/target/release/smol"
EDIT = f"{H}/smol-src/tools/ckpt_edit.py"
BASE = f"{H}/smol/restore/pub/box/seed1101/run.conf"
os.makedirs(out, exist_ok=True)

conf = []
for l in open(BASE):
    k = l.split()[0] if l.split() else ""
    if k in ("RANDSEED", "NSTEPS", "OFFSPRINGLOG", "SNAPEVERY", "CKPTEVERY", "KEEPCKPT", "CKPTKEEPEVERY"):
        continue
    conf.append(l.rstrip("\n"))
conf = [f"RANDSEED {seed}", f"NSTEPS {K}", "OFFSPRINGLOG 0", "SNAPEVERY 100000000", f"CKPTEVERY {K}",
        "KEEPCKPT 2", "CKPTKEEPEVERY 100000", "SPECIESFORGET 1", "COPYEVERY 1000", "COMPAT_DELETION 1"] + conf
open(f"{out}/run.conf", "w").write("\n".join(conf) + "\n")

def run(args):
    with open(f"{out}/log.txt", "a") as lg:
        return subprocess.run(args, cwd=out, stdout=lg, stderr=subprocess.STDOUT).returncode

copier, copied, born = set(), set(), {}
copies_off = 0; species_off = 0
def update_logs():
    global copies_off, species_off
    with open(f"{out}/copies.tsv") as f:
        f.seek(copies_off)
        for l in f:
            if not l.endswith("\n"): break
            copies_off += len(l)
            p = l.split("\t")
            if len(p) < 5 or not p[0].isdigit(): continue
            a, pa, c = int(p[1]), int(p[2]), int(p[3])
            if c == pa or c == a:
                copier.add(a); copied.add(c)
    with open(f"{out}/species.tsv") as f:
        f.seek(species_off)
        for l in f:
            if not l.endswith("\n"): break
            species_off += len(l)
            p = l.split("\t")
            if p[0].isdigit():
                born.setdefault(int(p[0]), int(p[1]))

def present(step):
    c = collections.Counter()
    with gzip.open(f"{out}/checkpoint_{step}.txt.gz", "rt") as f:
        for l in f:
            if l.startswith("MOLS"): break
        for l in f:
            c[int(l.split()[3])] += 1
    return c

tsv = open(f"{out}/suppress.tsv", "a")
if run([SMOL, "run", "run.conf"]) != 0:
    sys.exit("first chunk failed")
t = K
while t < nsteps:
    ck = f"checkpoint_{t}.txt.gz"
    if not os.path.exists(f"{out}/{ck}"):
        tsv.write(f"{t}\textinct\n"); tsv.flush(); break
    update_logs()
    mols = present(t); tot = sum(mols.values())
    if tot == 0:
        tsv.write(f"{t}\t0\t0\t0\n"); tsv.flush(); break
    P = [s for s in mols if s in copied and s not in copier and t - born.get(s, t) >= K]
    nP = sum(mols[s] for s in P)
    ops = []
    if mode == "P":
        for s in P: ops += ["--remove-species", str(s)]
    elif mode == "cull" and nP:
        ops = ["--cull", f"{1 - nP / tot:.6f}"]
    ops += ["--set", "NSTEPS", str(min(t + K, nsteps)), "--live-species", "--rng-seed", str(seed * 100000 + t // K)]
    r = subprocess.run(["python3", EDIT, ".", ck, "."] + ops, cwd=out, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit("edit failed: " + r.stderr)
    after = sum(present(t).values())
    tsv.write(f"{t}\t{tot}\t{nP}\t{tot - after}\n"); tsv.flush()
    # species.tsv was rewritten (live species only): restart the offset from its new end
    species_off = os.path.getsize(f"{out}/species.tsv")
    if run([SMOL, "resume", ck]) != 0:
        sys.exit("resume failed")
    for old in glob.glob(f"{out}/checkpoint_*.txt.gz"):
        s = int(old.rsplit("_", 1)[1].split(".")[0])
        if s <= t and s % 100000 != 0:
            os.remove(old)
    t += K
open(f"{out}/done", "w").write("ok\n")
