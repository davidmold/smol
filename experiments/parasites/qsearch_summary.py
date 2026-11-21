#!/usr/bin/env python3
"""qsearch_summary.py SET... -> per run set (e.g. ctl, qsearch), every 200k from the checkpoints:
runs alive, population, share of molecules carrying a label test (`?` + 2 or more letters), and
share of runs whose commonest species carries one. Parasite share from the time-shift windows
(OUT/sN/wNk/classes.tsv) when given as SET:OUT."""
import sys, os, re, gzip, glob, collections, statistics as S
LAB = re.compile(r"\?[A-Z]{2,}")
def mols(ck):
    c = collections.Counter()
    with gzip.open(ck, "rt") as f:
        for l in f:
            if l.startswith("MOLS"): break
        for l in f:
            p = l.split()
            if len(p) > 4: c[p[4].replace(".", "")] += 1
    return c
print("set        step   alive  pop(med)  label-test share of molecules (med)  top species has one  parasite share (med)")
for arg in sys.argv[1:]:
    st, _, ts = arg.partition(":")
    for step in range(200000, 2000001, 200000):
        pops, shares, tops, par = [], [], 0, []
        runs = sorted(glob.glob(f"{st}/s*"))
        for r in runs:
            ck = f"{r}/checkpoint_{step}.txt.gz"
            if not os.path.exists(ck): continue
            c = mols(ck); n = sum(c.values())
            if not n: continue
            pops.append(n); shares.append(sum(v for s, v in c.items() if LAB.search(s)) / n)
            tops += bool(LAB.search(c.most_common(1)[0][0]))
            if ts:
                cf = f"{ts}/{os.path.basename(r)}/w{step // 1000}k/classes.tsv"
                if os.path.exists(cf):
                    cl = [l.split() for l in open(cf)]; tot = sum(int(x[2]) for x in cl)
                    if tot: par.append(sum(int(x[2]) for x in cl if x[1] == "P") / tot)
        if pops:
            print(f"{os.path.basename(st):8} {step:>8}  {len(pops):>2}/{len(runs)}  {S.median(pops):7.0f}   {S.median(shares):.3f}"
                  f"                                  {tops}/{len(pops)}              {S.median(par) if par else float('nan'):.2f}")
