#!/usr/bin/env python3
"""Functional assay of a population: react its most abundant species with each other in the
reaction probe and summarise what the community does.

    assay.py RUN_DIR [--steps 100000,200000,...|--every 100000] [--top 12] [--strand]
             [--react PATH_TO_smol-react] [--max 2000000]

For each sampled step: the top-K species (by count) are reacted pairwise, each as initiator,
mutation off (smol-react). For every unordered pair {a,b} (a may equal b) with binding
probability > 0, a species counts as copied if either reaction produced an exact copy of it.
Pairs are weighted by abundance_a * abundance_b * P(bind). Reported per step:
  cover     share of molecules in the top-K
  selfrep   weighted share of same-species pairs that copy themselves (auto-replication)
  mutual    weighted share of productive cross pairs where each copies the other
  parasitic weighted share of productive cross pairs where only one is copied
            (Hickinbotham et al. 2021: R replicates P but P does not replicate R)
  rtime     mean probe steps of productive reactions (weighted)
  toplen, toprt  length of the most abundant species and its self-reaction steps
"""
import argparse, gzip, os, re, subprocess, tempfile

ap = argparse.ArgumentParser()
ap.add_argument("run"); ap.add_argument("--steps"); ap.add_argument("--every", type=int, default=100000)
ap.add_argument("--top", type=int, default=12); ap.add_argument("--strand", action="store_true")
ap.add_argument("--max", type=int, default=10**9)
ap.add_argument("--react", default=os.path.expanduser("~/smol/target/release/smol-react"))
a = ap.parse_args()

def ropen(p):
    return open(p) if os.path.exists(p) else gzip.open(p + ".gz", "rt")

seq = {}
for l in ropen(os.path.join(a.run, "species.tsv")):
    p = l.rstrip("\n").split("\t")
    if p[0].isdigit():
        seq[int(p[0])] = p[-1]
want = set(int(s) for s in a.steps.split(",")) if a.steps else None
pops = {}
cur = None
for l in ropen(os.path.join(a.run, "popdy001.dat")):
    t, s, c = l.split(","); t = int(t)
    if t > a.max:
        break
    if (want and t in want) or (not want and t % a.every == 0 and t > 0):
        pops.setdefault(t, {})[int(s)] = pops.get(t, {}).get(int(s), 0) + int(c)

pat = re.compile(r"(\S+)\s+initiates with (\S+)\s+P\(bind\)=([\d.]+)\s+copier=\S+\s+copied=\S+\s+steps=(\d+)")
print(f"== {a.run}")
print(f"{'step':>8} {'N':>6} {'cover':>5} {'selfrep':>7} {'mutual':>6} {'parasitic':>9} {'rtime':>6} {'toplen':>6} {'toprt':>6}")
for t in sorted(pops):
    pop = pops[t]; n = sum(pop.values())
    top = sorted(pop.items(), key=lambda kv: -kv[1])[:a.top]
    names = {f"s{sp}": sp for sp, _ in top}
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        for nm, sp in names.items():
            f.write(f"R {nm} {seq[sp]}\n")
        fn = f.name
    cmd = [a.react, fn, "5000"] + (["--strand"] if a.strand else [])
    out = subprocess.run(cmd, capture_output=True, text=True).stdout.split("\n")
    os.unlink(fn)
    res = {}; cur = None
    for l in out:
        m = pat.match(l)
        if m:
            cur = (m.group(1), m.group(2)); res[cur] = [float(m.group(3)), int(m.group(4)), set()]
        elif "-> product" in l and cur:
            res[cur][2].add(l.split()[2])
    w = lambda x: pop[names[x]] / n
    selfw = selfrep = mutual = para = prod_w = rt = rtw = 0.0
    keys = list(names)
    for i, x in enumerate(keys):
        for y in keys[i:]:
            r1 = res.get((x, y)); r2 = res.get((y, x))
            if not r1 or not r2:
                continue
            pb = max(r1[0], r2[0])
            if pb <= 0:
                continue
            wt = w(x) * w(y) * pb
            made = r1[2] | r2[2]
            productive = bool(made & {x, y})
            if productive:
                rt += wt * (r1[1] + r2[1]) / 2; rtw += wt
            if x == y:
                selfw += wt; selfrep += wt * (x in made)
            elif productive:
                prod_w += wt
                if x in made and y in made:
                    mutual += wt
                else:
                    para += wt
    t0 = keys[0]
    tr = res.get((t0, t0))
    toprt = tr[1] if tr and (t0 in tr[2]) else -1
    f2 = lambda v, d: v / d if d else float("nan")
    print(f"{t:>8} {n:>6} {sum(c for _, c in top)/n:>5.2f} {f2(selfrep, selfw):>7.2f} {f2(mutual, prod_w):>6.2f} {f2(para, prod_w):>9.2f} {f2(rt, rtw):>6.0f} {len(seq[names[t0]]):>6} {toprt:>6}")
