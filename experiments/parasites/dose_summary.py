#!/usr/bin/env python3
"""dose_summary.py DIR: summarise the graded PARASITEFREE 2 runs (dose.sh outputs in DIR).
q = probability that an exact copy of a non-self-copier is blocked: q000 = controls, q100 = pf2."""
import sys, re, collections, statistics as S
D = sys.argv[1]
levels = ["q000", "q025", "q050", "q075", "q100"]

surv = collections.defaultdict(list)
for l in open(f"{D}/dose_survival.txt"):
    d, rest = l.split(" ", 1)
    q = d.split("/")[1]
    surv[q].append("EXTINCT" not in rest)

first = collections.defaultdict(list)
cur = None
for l in open(f"{D}/dose_scanpath.txt"):
    if l.startswith("==="): cur = l.split()[1].split("/")[1]; continue
    m = re.search(r"first self-scanner: species \d+ at (\d+)", l)
    if m: first[cur].append(int(m.group(1)))
    elif "no self-scanner" in l: first[cur].append(None)

scan = collections.defaultdict(lambda: collections.defaultdict(list))
for l in open(f"{D}/dose_expo.txt"):
    p = l.split()
    if len(p) >= 5 and p[0] in levels:
        scan[p[0]][int(p[2])].append(float(p[4]) >= 1.5)

comp = {}
sect = False
for l in open(f"{D}/dose_complexity.txt"):
    if l.startswith("== medians"): sect = True; continue
    p = l.split()
    if sect and len(p) == 9 and p[0] in levels:
        comp[(p[0], int(p[1]))] = p[2:]

load = collections.defaultdict(list)
for l in open(f"{D}/dose_load.txt"):
    f = l.split() + ["0", "0"]
    d, b, ps = f[0], f[1], f[2]
    q = d.split("/")[1]
    b, ps = float(b), float(ps)
    if b + ps: load[q].append(ps / (b + ps))

print("q     survive  first self-scan (median; runs with one)   top replicator self-scans: share of snapshots 0.2-2M  [at 0.4M, 1M, 2M]")
for q in levels:
    f = [x for x in first[q] if x is not None]
    allsnap = [v for st in scan[q] for v in scan[q][st]]
    at = lambda st: f"{sum(scan[q][st])}/{len(scan[q][st])}" if scan[q][st] else "-"
    print(f"{q}  {sum(surv[q])}/{len(surv[q])}      {S.median(f) if f else float('nan'):>9.0f} ({len(f)}/{len(first[q])})"
          f"                         {sum(allsnap) / len(allsnap) if allsnap else float('nan'):.2f} ({sum(allsnap)}/{len(allsnap)})  [{at(400000)}, {at(1000000)}, {at(2000000)}]")
print("\nq     exact copies of non-self-copiers let through (share, median per run)")
for q in levels[1:]:
    print(f"{q}  {S.median(load[q]):.2f}" if load[q] else f"{q}  -")
print("\nComplexity medians (top replicator): runs, pop, len, '?', ops, self-reaction steps")
for st in (200000, 600000, 1000000, 2000000):
    for q in levels:
        if (q, st) in comp:
            n, pop, ln, qq, ops, rt, sel = comp[(q, st)]
            print(f"  {st:>8} {q}  n={n:>2} pop {pop:>6} len {ln:>3} ? {qq:>4} ops {ops:>5} rt {rt:>4}")
