#!/usr/bin/env python3
"""gm_summary.py SET TS -> gene-and-machine runs (SET/sN) every 200k: runs alive, population,
class shares from the time-shift windows (TS/sN/wNk/classes.tsv: R machines, G working genes,
P gene-parasites, N other), and the commonest machine's length, number of '?' and whether it
carries a label test ('?' + 2 or more letters)."""
import sys, os, re, glob, statistics as S, collections
SET, TS = sys.argv[1], sys.argv[2]
LAB = re.compile(r"\?[A-Z]{2,}")
print("step      alive  pop(med)   R     G     P     N    machine len  '?'  label-test  (medians; top machine of run s301)")
for step in range(200000, 2000001, 200000):
    sh = collections.defaultdict(list); pops = []; lens = []; qs = []; labs = 0; n = 0; ex = ""
    for r in sorted(glob.glob(f"{SET}/s*")):
        w = f"{TS}/{os.path.basename(r)}/w{step // 1000}k"
        if not os.path.exists(f"{w}/classes.tsv"): continue
        seq = {}
        for l in open(f"{w}/species.tsv") if os.path.exists(f"{w}/species.tsv") else []:
            p = l.rstrip("\n").split("\t")
            if p[0].isdigit() and len(p) >= 5: seq.setdefault(int(p[0]), p[-1])
        cl = [(int(a), b, int(c)) for a, b, c in (l.split() for l in open(f"{w}/classes.tsv"))]
        tot = sum(c for _, _, c in cl)
        if not tot: continue
        n += 1; pops.append(tot)
        for k in "RGPN": sh[k].append(sum(c for _, b, c in cl if b == k) / tot)
        reps = [(c, a) for a, b, c in cl if b == "R" and a in seq]
        if reps:
            m = seq[max(reps)[1]]
            lens.append(len(m)); qs.append(m.count("?")); labs += bool(LAB.search(m))
            if os.path.basename(r) == "s301": ex = m
    if n:
        print(f"{step:>8}  {n:>2}  {S.median(pops):7.0f}  " + "  ".join(f"{S.median(sh[k]):.2f}" for k in "RGPN") +
              f"   {S.median(lens) if lens else 0:5.0f}     {S.median(qs) if qs else 0:3.0f}    {labs}/{len(lens)}   {ex}")
