#!/usr/bin/env python3
"""Compare the dominant replicators of many runs.

    replicators.py RUN_DIR... [--window 200000] [--top 3]

For each run: over the last --window steps, a species is a replicator if, in copies.tsv,
it made exact copies (as copier) and was exactly copied at least 3 times each. Prints the
--top most abundant replicators of each run (abundance summed over popdy samples every
10k steps), then clusters every run's top replicator across runs by edit distance.
"""
import argparse, collections, os, sys
from rapidfuzz.distance import Levenshtein

ap = argparse.ArgumentParser()
ap.add_argument("runs", nargs="+")
ap.add_argument("--window", type=int, default=200000)
ap.add_argument("--top", type=int, default=3)
ap.add_argument("--thresh", type=float, default=0.3)
a = ap.parse_args()

tops = []
for run in a.runs:
    name = "/".join(run.rstrip("/").split("/")[-2:])
    seq = {}
    for l in open(os.path.join(run, "species.tsv")):
        p = l.rstrip("\n").split("\t")
        if len(p) >= 5:
            seq[int(p[0])] = p[4]
    last = 0
    abund = collections.Counter()
    rows = []
    with open(os.path.join(run, "popdy001.dat")) as f:
        for l in f:
            t, s, n = l.split(",")
            t = int(t)
            last = max(last, t)
            if t % 10000 == 0:
                rows.append((t, int(s), int(n)))
    lo = last - a.window
    for t, s, n in rows:
        if t > lo:
            abund[s] += n
    copier = collections.Counter(); copied = collections.Counter()
    for l in open(os.path.join(run, "copies.tsv")):
        p = l.split("\t")
        if int(p[0]) <= lo:
            continue
        act, pas, prod, n = int(p[1]), int(p[2]), int(p[3]), int(p[4])
        if prod == pas:
            copier[act] += n; copied[pas] += n
        elif prod == act:
            copier[act] += n; copied[act] += n
    reps = [s for s in abund if copier[s] >= 3 and copied[s] >= 3]
    reps.sort(key=lambda s: -abund[s])
    tot = sum(abund.values()) or 1
    print(f"== {name}  (last step {last})")
    for s in reps[: a.top]:
        print(f"   {100 * abund[s] / tot:5.1f}%  len {len(seq[s]):3d}  {seq[s]}")
    if reps:
        tops.append((name, seq[reps[0]], 100 * abund[reps[0]] / tot))

print(f"\nTop replicator of each run, clustered (edit distance <= {a.thresh}):")
n = len(tops)
parent = list(range(n))
def find(i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]; i = parent[i]
    return i
for i in range(n):
    for j in range(i + 1, n):
        if Levenshtein.normalized_distance(tops[i][1], tops[j][1]) <= a.thresh:
            parent[find(i)] = find(j)
groups = collections.defaultdict(list)
for i in range(n):
    groups[find(i)].append(i)
for g in sorted(groups.values(), key=len, reverse=True):
    print(f"\n  cluster of {len(g)}:")
    for i in g:
        print(f"    {tops[i][0]:<22} {tops[i][2]:5.1f}%  {tops[i][1]}")
