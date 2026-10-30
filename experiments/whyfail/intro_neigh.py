#!/usr/bin/env python3
"""Test on existing introductions: were the W26 (and S21, E1) introductions that established
surrounded by more relatives? (run on aida)

For each natural.py run into seed21021 at 9.9M / 10.0M: the introduced cells (intro.txt), and,
in the recipient checkpoint (runs/natural/_src/<rec>), the molecules in the 8 cells around them
(and in the 5x5 block). A neighbour is a 'relative' if its sequence carries the S21-family copy
core '>C$==' (W26: Q>BO^$O%H>C$===?O>A$$T$TQ>). Outcome = established (endpoint >= 10%,
introduced family) from experiments/processive/natural_families.txt.
"""
import gzip, json, os, re, statistics as st
from collections import defaultdict

H = os.path.expanduser("~/smol")
fam = {}
for l in open(f"{H}/experiments/processive/natural_families.txt"):
    p = l.rstrip("\n").split("\t")
    if p[0] != "run" and p[3] == "introduced":
        fam[p[0]] = p[2]
def grid(rec):
    step = {"s21_9.9M": 9900000, "s21_10.0M": 10000000}[rec]
    L = gzip.open(f"{H}/runs/natural/_src/{rec}/checkpoint_{step}.txt.gz", "rt").read().split("\n")
    n = int(L[4].split()[1]); conf = L[5:5 + n]
    gx = next(int(float(l.split()[1])) for l in conf if l.split()[:1] == ["GRIDX"])
    gy = next(int(float(l.split()[1])) for l in conf if l.split()[:1] == ["GRIDY"])
    i = 5 + n; m = int(L[i + 2].split()[1]); g = {}
    for l in L[i + 3:i + 3 + m]:
        p = l.split(); g[(int(p[1]), int(p[2]))] = p[4].replace(".", "")
    return gx, gy, g
rel = lambda s: ">C$==" in s
out = defaultdict(list)
for rec in ("s21_9.9M", "s21_10.0M"):
    gx, gy, g = grid(rec)
    for seq in ("W26", "S21", "E1"):
        for n in (1, 2):
            for k in range(16):
                d = f"{H}/runs/natural/{rec}/{seq}/n{n}/rep{k}"
                if not os.path.exists(f"{d}/intro.txt"):
                    continue
                cells = [tuple(map(int, c)) for c in re.findall(r"\((\d+), (\d+)\)", open(f"{d}/intro.txt").read())]
                nb8 = set(); nb24 = set()
                for (x, y) in cells:
                    for dx in range(-2, 3):
                        for dy in range(-2, 3):
                            c = ((x + dx) % gx, (y + dy) % gy)
                            if c in cells:
                                continue
                            nb24.add(c)
                            if abs(dx) <= 1 and abs(dy) <= 1:
                                nb8.add(c)
                r8 = sum(1 for c in nb8 if c in g and rel(g[c])); r24 = sum(1 for c in nb24 if c in g and rel(g[c]))
                e8 = sum(1 for c in nb8 if c not in g)
                key = f"{rec}/{seq}/n{n}/rep{k}"
                out[(rec, seq)].append(dict(run=key, est=key in fam, rel8=r8, rel24=r24, empty8=e8, n8=len(nb8)))
print(f"{'recipient':10} {'seq':4} {'est':>4} {'rel8 est':>9} {'rel8 fail':>10} {'rel24 est':>10} {'rel24 fail':>11} {'est | rel8>=2':>14} {'est | rel8<2':>13}")
allrows = []
for (rec, seq), v in sorted(out.items()):
    E = [r for r in v if r["est"]]; F = [r for r in v if not r["est"]]; allrows += v
    m = lambda xs, f: f"{st.mean(f(r) for r in xs):.2f}" if xs else "-"
    hi = [r for r in v if r["rel8"] >= 2]; lo = [r for r in v if r["rel8"] < 2]
    print(f"{rec:10} {seq:4} {len(E):>4} {m(E, lambda r: r['rel8']):>9} {m(F, lambda r: r['rel8']):>10} "
          f"{m(E, lambda r: r['rel24']):>10} {m(F, lambda r: r['rel24']):>11} "
          f"{sum(r['est'] for r in hi):>6}/{len(hi):<7} {sum(r['est'] for r in lo):>5}/{len(lo):<7}")
json.dump(allrows, open("intro_neigh.json", "w"), indent=1)
# the recipient's overall relative density, for reference
for rec in ("s21_9.9M", "s21_10.0M"):
    gx, gy, g = grid(rec)
    print(f"{rec}: {len(g)} molecules, relatives (carry '>C$==') {sum(1 for s in g.values() if rel(s))} "
          f"({sum(1 for s in g.values() if rel(s)) / len(g):.2%})")
