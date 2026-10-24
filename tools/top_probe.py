#!/usr/bin/env python3
"""Top species at chosen steps, each reacted with itself in the probe (classic chemistry).

    top_probe.py RUN_DIR STEP[,STEP...] [--top 6] [--strand]

For each species: share, length, and in a self-reaction: steps (5000 = still running at the
limit), number of exact self-copies made, and other products. A processive replicator makes
many copies in one reaction.
"""
import argparse, gzip, os, re, subprocess, tempfile

ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("steps")
ap.add_argument("--top", type=int, default=6); ap.add_argument("--strand", action="store_true")
ap.add_argument("--react", default=os.path.expanduser("~/smol/target/release/smol-react"))
a = ap.parse_args()
ropen = lambda p: open(p) if os.path.exists(p) else gzip.open(p + ".gz", "rt")
want = sorted(int(s) for s in a.steps.split(",")); W = set(want)
pops = {}
for l in ropen(f"{a.run}/popdy001.dat"):
    t = int(l[:l.index(",")])
    if t in W:
        _, s, c = l.split(","); pops.setdefault(t, {})[int(s)] = int(c)
    elif t > want[-1]:
        break
need = {s for p in pops.values() for s, _ in sorted(p.items(), key=lambda kv: -kv[1])[:a.top]}
seq = {}
for l in ropen(f"{a.run}/species.tsv"):
    i = l.find("\t")
    if i > 0 and l[:i].isdigit() and int(l[:i]) in need:
        p = l.rstrip("\n").split("\t"); seq[int(p[0])] = (p[-1], int(p[1]))
for t in want:
    pop = pops.get(t, {}); n = sum(pop.values())
    print(f"== step {t}: {n} molecules")
    for s, c in sorted(pop.items(), key=lambda kv: -kv[1])[:a.top]:
        sq, born = seq[s]
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
            f.write(f"R x {sq}\n"); fn = f.name
        out = subprocess.run([a.react, fn, "5000"] + (["--strand"] if a.strand else []), capture_output=True, text=True).stdout
        os.unlink(fn)
        m = re.search(r"steps=(\d+)", out); steps = int(m.group(1)) if m else -1
        prods = re.findall(r"-> product (\S+)", out)
        print(f"  {c/n:.3f} born {born:>8} len {len(sq):>3} self-reaction steps {steps:>4} copies {prods.count('x'):>3} other {len(prods)-prods.count('x'):>2}  {sq[:60]}")
