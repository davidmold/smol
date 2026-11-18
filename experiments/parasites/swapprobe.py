#!/usr/bin/env python3
"""swapprobe.py -> does the 2M top replicator handle the 2M parasites differently from the 1M one?

For each published run: top replicator at 1M (R1) and at 2M (R2) (decoy.py's replicator
definition), and the 2M parasites. Each replicator is reacted with each parasite (smol-react,
--compat); bound time (parasite abundance x P(bind) x steps) is split into copied / held >= 50
steps with no copy / brief. If R2 holds where R1 copies, that explains why 1M replicators put in
2M worlds keep a smaller population (findings 20, swap experiment)."""
import os, collections
g = {"__file__": os.path.join(os.path.dirname(os.path.abspath(__file__)), "decoy.py")}
src = open(g["__file__"]).read().split('if mode == "survey":')[0]
src = src.replace("mode, OUT = sys.argv[1], os.path.abspath(sys.argv[2])", "mode, OUT = 'lib', '/tmp'")
exec(src, g)
C, analyse, react_all = g["C"], g["analyse"], g["react_all"]

def reps(c):
    top, ab, selfc, paras, t, dtrap, decoys = analyse(c)
    return top, ab, [i for i in range(len(top)) if selfc[i]], paras

print("run   host     vs 2M parasites:  copied  held  brief")
summ = collections.defaultdict(list)
for r in range(1, 9):
    t1, a1, r1, _ = reps(C.snapshot(r, 1000000))
    t2, a2, r2, p2 = reps(C.snapshot(r, 2000000))
    if not r1 or not r2 or not p2:
        print(f"run{r}  skipped (no replicator or parasite found)"); continue
    R1 = max(r1, key=lambda i: a1[i]); R2 = max(r2, key=lambda i: a2[i])
    for lab, R in (("1M rep", t1[R1]), ("2M rep", t2[R2])):
        res = react_all([R] + [t2[p] for p in p2])
        tc = collections.Counter()
        for k, p in enumerate(p2, 1):
            pb, st, prods = res.get((0, k), [0, 0, set()])[:3]
            w = a2[p] * pb * st
            tc["copied" if k in prods else ("held" if st >= 50 else "brief")] += w
        tot = sum(tc.values()) or 1
        summ[lab].append((tc["copied"] / tot, tc["held"] / tot))
        print(f"run{r}  {lab}                    {tc['copied'] / tot:5.2f}   {tc['held'] / tot:5.2f}  {tc['brief'] / tot:5.2f}")
import statistics as S
for lab, v in summ.items():
    print(f"median {lab}: copied {S.median(x for x, _ in v):.2f}, held {S.median(y for _, y in v):.2f}")
