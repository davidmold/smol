#!/usr/bin/env python3
"""agg.py DIR STEP: population after each treatment (mean of 3 seeds) per run, at several
times after the edit, from log.txt ("Step N done, number of molecules is X")."""
import sys, re, glob, os, statistics, collections
D, S = sys.argv[1], int(sys.argv[2])
offs = [0, 2000, 10000, 20000, 50000, 100000]
T = ["base", "noP", "cullP", "noN", "cullN"]
pop = collections.defaultdict(list)   # (run, t, off) -> [pop]
for f in glob.glob(f"{D}/r*/*/s*/log.txt"):
    r, t = f.split("/")[-4], f.split("/")[-3]
    vals = {}
    for l in open(f):
        m = re.match(r"Step (\d+) done, number of molecules is (\d+)", l)
        if m: vals[int(m.group(1))] = int(m.group(2))
    for o in offs:
        if S + o in vals: pop[(r, t, o)].append(vals[S + o])
runs = sorted({k[0] for k in pop}, key=lambda x: int(x[1:]))
for o in offs:
    print(f"== {o} steps after the edit (mean of seeds; molecules)")
    print("run    " + "".join(f"{t:>8}" for t in T) + "   noP-cullP  noN-cullN")
    dP, dN = [], []
    for r in runs:
        m = {t: (statistics.mean(pop[(r, t, o)]) if pop[(r, t, o)] else float('nan')) for t in T}
        dP.append(m["noP"] - m["cullP"]); dN.append(m["noN"] - m["cullN"])
        print(f"{r:6} " + "".join(f"{m[t]:8.0f}" for t in T) + f"   {dP[-1]:9.0f}  {dN[-1]:9.0f}")
    ok = [x for x in dP if x == x]; okn = [x for x in dN if x == x]
    if ok: print(f"median  noP-cullP {statistics.median(ok):.0f} (positive in {sum(x>0 for x in ok)}/{len(ok)} runs);  noN-cullN {statistics.median(okn):.0f} (positive in {sum(x>0 for x in okn)}/{len(okn)})")
