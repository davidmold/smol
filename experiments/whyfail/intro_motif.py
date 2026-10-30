#!/usr/bin/env python3
"""Introduction outcome vs relatives next to the introduced cells, for nested definitions of
'relative' (run on aida in ~/published/whyfail after intro_neigh.py)."""
import gzip, json, math, os, re

H = os.path.expanduser("~/smol")
rows = json.load(open("intro_neigh.json"))
def grid(rec):
    step = {"s21_9.9M": 9900000, "s21_10.0M": 10000000}[rec]
    L = gzip.open(f"{H}/runs/natural/_src/{rec}/checkpoint_{step}.txt.gz", "rt").read().split("\n")
    n = int(L[4].split()[1]); i = 5 + n; m = int(L[i + 2].split()[1]); g = {}
    for l in L[i + 3:i + 3 + m]:
        p = l.split(); g[(int(p[1]), int(p[2]))] = p[4].replace(".", "")
    return g
G = {r: grid(r) for r in ("s21_9.9M", "s21_10.0M")}

def mw(a, b):
    """One-sided Mann-Whitney (a larger), normal approximation with tie-averaged ranks."""
    n1, n2 = len(a), len(b); allv = sorted(a + b); ranks = {}
    for v in set(allv):
        idx = [i + 1 for i, x in enumerate(allv) if x == v]; ranks[v] = sum(idx) / len(idx)
    U = sum(ranks[x] for x in a) - n1 * (n1 + 1) / 2
    mu = n1 * n2 / 2; sd = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    return U, (0.5 * math.erfc(((U - mu - 0.5) / sd) / math.sqrt(2)) if sd else 1.0)

motifs = [("W26 motif 'O%H>C$===?O>A'", "O%H>C$===?O>A"), ("W26 core '>C$===?O>A'", ">C$===?O>A"),
          ("family '>C$=='", ">C$==")]
import sys
ONLY_N = sys.argv[1] if len(sys.argv) > 1 else None   # e.g. n2
print("fraction of neighbour cells holding a relative" + (f", {ONLY_N} only" if ONLY_N else ""))
for name, mo in motifs:
    g0 = G["s21_9.9M"]
    print(f"{name}: carried by {sum(mo in s for s in g0.values()) / len(g0):.1%} of molecules at 9.9M")
    for seq in ("W26", "S21", "E1"):
        E, F = [], []
        for r in rows:
            rec, s, n, rep = r["run"].split("/")
            if s != seq:
                continue
            cells = [(int(a), int(b)) for a, b in re.findall(r"\((\d+), (\d+)\)", open(f"{H}/runs/natural/{r['run']}/intro.txt").read())]
            nb = {((x + dx) % 125, (y + dy) % 100) for x, y in cells for dx in (-1, 0, 1) for dy in (-1, 0, 1)} - set(cells)
            c = sum(1 for q in nb if q in G[rec] and mo in G[rec][q]) / len(nb)  # fraction of neighbour cells
            if ONLY_N and n != ONLY_N:
                continue
            (E if r["est"] else F).append(round(c, 2))
        if E and F:
            U, p = mw(E, F)
            print(f"   {seq}: established {sorted(E)}  failed mean {sum(F) / len(F):.2f} (n={len(F)})  one-sided MW p = {p:.3f}")
        else:
            print(f"   {seq}: established {len(E)}; failed mean {(sum(F) / len(F)) if F else 0:.2f} (n={len(F)})")
