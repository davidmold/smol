#!/usr/bin/env python3
"""swap_agg.py DIR: the 2x2 time-shift of replicators (R) and parasites (P) between each run's
1M and 2M worlds. Every edit replaces one class with molecules drawn from the same run at the
same (control) or the other time. Per run, mean of 5 seeds, at several times after the edit:
population, and parasite share = molecules whose sequence was classed P at 1M or 2M.

  1M world:  self1 / pself1 = R1+P1 (controls)   late2 = R2+P1   plate2 = R1+P2
  2M world:  self2 / pself2 = R2+P2 (controls)   early1 = R1+P2  pearly1 = R2+P1
"""
import sys, re, glob, statistics, collections, os
D = sys.argv[1]; P = os.path.dirname(D.rstrip("/"))
offs = [2000, 10000, 20000, 50000, 100000]
def pseqs(src):
    sq = {}
    for l in open(f"{src}/species.tsv"):
        p = l.rstrip("\n").split("\t")
        if p[0].isdigit(): sq.setdefault(int(p[0]), p[-1])
    return {sq[int(l.split()[0])] for l in open(f"{src}/classes.tsv") if l.split()[1] == "P" and int(l.split()[0]) in sq}
PS = {}
pop = collections.defaultdict(list); psh = collections.defaultdict(list)
for f in glob.glob(f"{D}/r*/*/s*/log.txt"):
    d = os.path.dirname(f); r, t = f.split("/")[-4], f.split("/")[-3]
    S = 1000000 if t in ("self1", "late2", "pself1", "plate2") else 2000000
    if r not in PS:
        PS[r] = pseqs(f"{P}/s1000000/{r}/src") | pseqs(f"{P}/s2000000/{r}/src")
    vals = {}
    for l in open(f):
        m = re.match(r"Step (\d+) done, number of molecules is (\d+)", l)
        if m: vals[int(m.group(1))] = int(m.group(2))
    if S + offs[-1] not in vals: continue      # unfinished
    sq = {}
    for l in open(f"{d}/species.tsv"):
        p = l.rstrip("\n").split("\t")
        if p[0].isdigit(): sq.setdefault(int(p[0]), p[-1])
    want = {S + o for o in offs}; tot = collections.Counter(); par = collections.Counter()
    for l in open(f"{d}/popdy001.dat"):
        a = l.split(","); s = int(a[0])
        if s in want:
            c = int(a[2]); tot[s] += c
            if sq.get(int(a[1])) in PS[r]: par[s] += c
    for o in offs:
        pop[(r, t, o)].append(vals[S + o])
        if tot[S + o]: psh[(r, t, o)].append(par[S + o] / tot[S + o])
runs = sorted({k[0] for k in pop}, key=lambda x: int(x[1:]))
mean = lambda v: statistics.mean(v) if v else float("nan")
comps = [("late2", "self1", "1M world: 2M replicators vs own"),
         ("plate2", "pself1", "1M world: 2M parasites vs own"),
         ("early1", "self2", "2M world: 1M replicators vs own"),
         ("pearly1", "pself2", "2M world: 1M parasites vs own")]
for a, b, label in comps:
    print(f"#### {label}   ({a} vs control {b})")
    for o in offs:
        d = []; dp = []
        rows = []
        for r in runs:
            x, y = mean(pop[(r, b, o)]), mean(pop[(r, a, o)]); px, py = mean(psh[(r, b, o)]), mean(psh[(r, a, o)])
            d.append(y - x); dp.append(py - px); rows.append(f"{r}:{y-x:+.0f}")
        ok = [v for v in d if v == v]; okp = [v for v in dp if v == v]
        if ok:
            print(f"  +{o:6d}  pop diff median {statistics.median(ok):+6.0f} (up in {sum(v>0 for v in ok)}/{len(ok)})"
                  f"   parasite share diff median {statistics.median(okp):+.2f} (up in {sum(v>0 for v in okp)}/{len(okp)})   " + " ".join(rows))
print("\nControls (mean population / parasite share at +10000):")
for r in runs:
    print("  " + r + "  " + "  ".join(f"{t}={mean(pop[(r,t,10000)]):.0f}/{mean(psh[(r,t,10000)]):.2f}" for t in ("self1","pself1","self2","pself2")))
