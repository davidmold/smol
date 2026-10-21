#!/usr/bin/env python3
"""Population size and molecule length over time, pooled across runs, for comparing run sets.

    length_traj.py LABEL RUN_DIR... [--every 100000] [--max 2000000]

Each RUN_DIR needs popdy001.dat and species.tsv (plain or .gz). For each sampled step, per run:
N molecules, median and 90th-percentile length (molecule-weighted), length and share of the
most abundant species. Prints the median of each across runs alive at that step.
"""
import gzip, os, statistics, sys

args = [a for a in sys.argv[1:]]
every, mx = 100000, 2000000
if "--every" in args:
    i = args.index("--every"); every = int(args[i + 1]); del args[i:i + 2]
if "--max" in args:
    i = args.index("--max"); mx = int(args[i + 1]); del args[i:i + 2]
label, runs = args[0], args[1:]

def ropen(p):
    return open(p) if os.path.exists(p) else gzip.open(p + ".gz", "rt")

def wq(pairs, q):  # weighted quantile of (value, weight)
    pairs = sorted(pairs); tot = sum(w for _, w in pairs); acc = 0
    for v, w in pairs:
        acc += w
        if acc >= q * tot:
            return v
    return pairs[-1][0]

table = {}
for d in runs:
    seq = {}
    for l in ropen(os.path.join(d, "species.tsv")):
        p = l.rstrip("\n").split("\t")
        if p[0].isdigit():
            seq[int(p[0])] = p[-1]
    cur = {}; step = None
    def flush(step, cur):
        if step is not None and step % every == 0 and step <= mx and cur:
            pairs = [(len(seq.get(s, "")), c) for s, c in cur.items()]
            n = sum(c for _, c in pairs)
            top = max(cur.items(), key=lambda kv: kv[1])
            table.setdefault(step, []).append((n, wq(pairs, .5), wq(pairs, .9), len(seq.get(top[0], "")), top[1] / n))
    for l in ropen(os.path.join(d, "popdy001.dat")):
        t, s, c = l.split(",")
        t = int(t)
        if t != step:
            flush(step, cur); step = t; cur = {}
        cur[int(s)] = cur.get(int(s), 0) + int(c)
    flush(step, cur)

print(f"== {label}: {len(runs)} runs (medians across runs alive at each step)")
print(f"{'step':>8} {'runs':>4} {'N':>6} {'medlen':>6} {'p90len':>6} {'toplen':>6} {'topshare':>8}")
for t in sorted(table):
    r = table[t]
    m = lambda k: statistics.median(x[k] for x in r)
    print(f"{t:>8} {len(r):>4} {m(0):>6.0f} {m(1):>6.0f} {m(2):>6.0f} {m(3):>6.0f} {m(4):>8.2f}")
