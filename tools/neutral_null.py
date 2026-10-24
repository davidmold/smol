#!/usr/bin/env python3
"""Neutral drift null for the family metrics.

    neutral_null.py RUN_DIR OUT_DIR [--start 300000] [--end 1000000] [--seed 1]

Starts from the run's actual population at --start and continues it as a neutral Moran
process: each step about N*DECAY molecules die at random and are replaced by copies of
randomly chosen survivors, with no regard to sequence. A copy is a new species with
probability r, the run's own rate of new species per replacement over the 100k steps
before --start, by one substitution, insertion or deletion at a random position.

Writes popdy001.dat and species.tsv in OUT_DIR in the simulator's formats, so
families.py and window_novelty.py run on it unchanged. Same population size, same
turnover, same rate of new sequences: only selection is missing.
"""
import argparse, os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from families import load_species, read_popdy

ap = argparse.ArgumentParser()
ap.add_argument("run"); ap.add_argument("out")
ap.add_argument("--start", type=int, default=300000)
ap.add_argument("--end", type=int, default=1000000)
ap.add_argument("--decay", type=float, default=0.0005)
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--placed", action="store_true", help="calibrate r from new species whose first product was placed (needs OFFSPRINGLOG copies.tsv)")
a = ap.parse_args()
rng = random.Random(a.seed)

seq, birth = load_species(a.run)
pop0 = read_popdy(a.run, 1000)
counts = pop0.get(a.start)
if not counts:
    sys.exit(f"no population at step {a.start}")
N = sum(counts.values())
if a.placed:
    import gzip as _gz
    _op = lambda p: open(p) if os.path.exists(p) else _gz.open(p + ".gz", "rt")
    first = {}
    for l in _op(os.path.join(a.run, "copies.tsv")):
        r = l.rstrip("\n").split("\t")
        if len(r) > 6 and int(r[6]) > 0:
            first.setdefault(int(r[3]), int(r[0]))
    # a species counts as born when its first placed product is logged (rows are per 1000 steps)
    born = sum(1 for s, t in first.items() if a.start - 100000 < t <= a.start and birth.get(s, 0) > a.start - 101000)
else:
    born = sum(1 for b in birth.values() if a.start - 100000 < b <= a.start)
r = born / (N * a.decay * 100000)
alphabet = sorted({ch for s in seq.values() for ch in s})
print(f"N={N} new species in previous 100k={born} r={r:.4f} alphabet={len(alphabet)}", flush=True)

pop = [s for s, c in counts.items() for _ in range(c)]
seq_id = {}
for s in set(pop):
    seq_id.setdefault(seq[s], s)
next_id = max(seq) + 1
os.makedirs(a.out, exist_ok=True)
sp_out = open(os.path.join(a.out, "species.tsv"), "w")
for s in sorted(set(pop)):
    sp_out.write(f"{s}\t{birth.get(s, 0)}\t-1\t-1\t{seq[s]}\n")
pd_out = open(os.path.join(a.out, "popdy001.dat"), "w")
# earlier history, so windows before --start look like the real run
for t, sp in pop0.items():
    if t < a.start:
        for s, c in sorted(sp.items()):
            pd_out.write(f"{t},{s},{c}\n")
            if s not in seq_id.values():
                pass
known = set(pop)
for t, sp in pop0.items():
    if t < a.start:
        for s in sp:
            if s not in known:
                known.add(s); sp_out.write(f"{s}\t{birth.get(s, 0)}\t-1\t-1\t{seq[s]}\n")

def mutate(x):
    k = rng.randrange(3); i = rng.randrange(len(x)) if x else 0
    if k == 0 and x:
        return x[:i] + rng.choice(alphabet) + x[i + 1:]
    if k == 1 or not x:
        return x[:i] + rng.choice(alphabet) + x[i:]
    return x[:i] + x[i + 1:] if len(x) > 1 else x

lam = N * a.decay
for t in range(a.start, a.end + 1):
    if t % 1000 == 0:
        c = {}
        for s in pop:
            c[s] = c.get(s, 0) + 1
        for s, n in sorted(c.items()):
            pd_out.write(f"{t},{s},{n}\n")
    # Poisson(lam) replacements; lam is ~5 so a simple loop is fine
    k = 0; p = rng.random(); L = 2.718281828 ** (-lam); acc = L; cum = L
    while p > cum:
        k += 1; acc *= lam / k; cum += acc
    for _ in range(k):
        d = rng.randrange(N); par = pop[rng.randrange(N)]
        if rng.random() < r:
            ns = mutate(seq[par])
            sid = seq_id.get(ns)
            if sid is None:
                sid = next_id; next_id += 1
                seq[sid] = ns; seq_id[ns] = sid
                sp_out.write(f"{sid}\t{t}\t{par}\t{par}\t{ns}\n")
            pop[d] = sid
        else:
            pop[d] = par
sp_out.close(); pd_out.close()
print("done", flush=True)
