#!/usr/bin/env python3
"""decoy.py survey OUT -> how common are decoys, across the paper's 8 runs and over time?

Decoy: a molecule that cannot copy itself but binds parasites and holds them without copying them
(run 4 of the competition, findings 26). For each run and each 200k snapshot of the published
runs: take the 30 commonest species; classify each as replicator (makes an exact copy of a partner or itself in some pair) or not;
react all pairs (smol-react, --compat). Parasites: species that never make a copy but are copied exactly by a
replicator. For every parasite p and partner X, the expected bound time per encounter is
P(bind) x steps, weighted by the abundances of p and X. That time is split into:
  copied   X makes an exact copy of p
  trapped  no copy of p, held >= TRAP steps; split by whether X is a self-copier (replicator
           trap) or not (decoy trap)
  brief    no copy, released within TRAP steps
Prints per snapshot: population, parasite share of the top 30, and the shares of parasite-bound
time copied / trapped by replicators / trapped by decoys / brief, plus the decoys' share of
molecules (non-self-copiers whose decoy-trap time exceeds the time they spend copying anything).
Also writes decoys.tsv (run, step, sequence, molecules, trap time) for the removal test.
"""
import sys, os, re, collections, types, subprocess, tempfile
H = os.path.expanduser("~")
C = types.SimpleNamespace()
_src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "compete.py")).read()
_src = _src.split('if mode == "setup":')[0].replace("mode, OUT = sys.argv[1], sys.argv[2]", "")
exec(_src, C.__dict__)
mode, OUT = sys.argv[1], os.path.abspath(sys.argv[2])
TRAP, K = 50, 30
RUNS = [int(x) for x in os.environ.get("RUNS", "1,2,3,4,5,6,7,8").split(",")]
STEPS = [int(x) for x in os.environ.get("STEPS", ",".join(str(s) for s in range(200000, 2000001, 200000))).split(",")]

def react_all(seqs):
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        for i, s in enumerate(seqs): f.write(f"R m{i} {s}\n")
        n = f.name
    out = subprocess.run([C.REACT, n, "5000", "--compat"], capture_output=True, text=True).stdout; os.unlink(n)
    res = {}; cur = None
    for l in out.splitlines():
        m = re.match(r"m(\d+)\s+initiates with m(\d+)\s+P\(bind\)=([\d.]+)\s+copier=m(\d+)\s+copied=m(\d+)\s+steps=(\d+)", l)
        if m:
            a, b = int(m.group(1)), int(m.group(2))
            cur = (min(a, b), max(a, b))
            if cur in res: cur = None; continue
            res[cur] = [float(m.group(3)), int(m.group(6)), set(), int(m.group(4)), int(m.group(5))]; continue
        m = re.match(r"\s+-> (?:product|passive made) m(\d+)", l)
        if m and cur: res[cur][2].add(int(m.group(1)))
    return res

def analyse(c):
    top = [s for s, _ in c.most_common(K)]
    ab = [c[s] for s in top]
    res = react_all(top)
    # replicator: makes an exact copy of its partner (or of itself) in some pair;
    # parasite: never does, but is copied exactly by a replicator (classify.py's R and P)
    n = len(top)
    selfc = [False] * n
    for (a, b), (pb, st, prods, cop, cpd) in res.items():
        if pb > 0 and cpd in prods: selfc[cop] = True
    paras = [i for i in range(n) if not selfc[i] and any(
        v[0] > 0 and v[4] == i and i in v[2] and selfc[v[3]] for v in res.values())]
    t = collections.Counter(); dtrap = collections.Counter(); dwork = collections.Counter()
    for p in paras:
        for x in range(len(top)):
            if x == p: continue
            pb, st, prods = res.get((min(p, x), max(p, x)), [0, 0, set()])[:3]
            w = ab[p] * ab[x] * pb * st
            if not w: continue
            if p in prods: t["copied"] += w
            elif st >= TRAP:
                if selfc[x]:
                    t["rtrap"] += w
                    t["rtrap_mut" if prods else "rtrap_none"] += w
                else: t["dtrap"] += w; dtrap[x] += w
            else: t["brief"] += w
    # decoys: non-self-copiers whose decoy-trap time exceeds the time they spend copying anything
    for x in range(len(top)):
        if selfc[x]: continue
        for y in range(len(top)):
            pb, st, prods = res.get((min(x, y), max(x, y)), [0, 0, set()])[:3]
            if y in prods and x != y: dwork[x] += ab[x] * ab[y] * pb * st
    decoys = [x for x in dtrap if dtrap[x] > dwork[x] and x not in paras]
    return top, ab, selfc, paras, t, dtrap, decoys

if mode == "survey":
    os.makedirs(OUT, exist_ok=True)
    dt = open(f"{OUT}/decoys.tsv", "a")
    print("run   step     pop  par%top   copied  rep-trap  decoy-trap  brief   rt:none rt:other-product   decoy-mol%  ndecoy", flush=True)
    for r in RUNS:
        for step in STEPS:
            try: c = C.snapshot(r, step)
            except FileNotFoundError: continue
            pop = sum(c.values())
            top, ab, selfc, paras, t, dtrap, decoys = analyse(c)
            tot = (t["copied"] + t["rtrap"] + t["dtrap"] + t["brief"]) or 1
            topn = sum(ab)
            print(f"run{r} {step:>8} {pop:>6}   {sum(ab[p] for p in paras) / topn:5.2f}    {t['copied'] / tot:5.2f}    {t['rtrap'] / tot:5.2f}      {t['dtrap'] / tot:5.2f}   {t['brief'] / tot:5.2f}   {t['rtrap_none'] / tot:5.2f} {t['rtrap_mut'] / tot:5.2f}   {sum(ab[x] for x in decoys) / pop * 100:5.1f}     {len(decoys):>3}", flush=True)
            for x in decoys:
                dt.write(f"run{r}\t{step}\t{top[x]}\t{ab[x]}\t{dtrap[x]:.0f}\n")
            dt.flush()
