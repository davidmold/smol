#!/usr/bin/env python3
"""Summarise partnership batch 2 (~/smol/runs/b2).

    agg2.py B2_DIR

invade/, reverse/, bottleneck/, rescue/: per condition, runs alive at the end and the median
final share of the introduced or tracked sequences (exact matches).
access/: per (seed, checkpoint time) and kind of run, how many runs had an established
partnership: in some 10k window after the start, >= 10% of placed exact copies happened in
mutual cross-species pairs (each direction >= 20 placed copies in that window). For intro runs
also the final count of the introduced sequences.
"""
import collections, glob, gzip, os, statistics, sys

B2 = sys.argv[1]
COP = "RERERERE^BQQ$=?>%}"
PAIRS = {"s1": ("S^BO$=J?>G}", "F^O>"), "s2": ("RER$Y$^BQQ$=?>H}", "ERE^AEBP>"),
         "s3": ("FSF^BQ$=?>G}", "SES^AQ"), "g1": ("FUS^YF^BPQ$=?L>}", "SIF^ABQ")}
OWN = {"18018": ("REER^BQQ$=?>G}", "ERTER^AQQ$T?>G}"), "20020": ("RECRERDRD^BQQ$=?>%}", "ERRE^AQQ"),
       "23023": ("EERRE^BPQ$=?>%}", "RRE^AQQ")}

def ropen(p):
    return open(p) if os.path.exists(p) else gzip.open(p + ".gz", "rt")

def species(d):
    s = {}
    for l in ropen(f"{d}/species.tsv"):
        p = l.rstrip("\n").split("\t")
        if p[0].isdigit():
            s[int(p[0])] = p[-1]
    return s

def final(d, tracked):
    if not (os.path.exists(f"{d}/species.tsv") or os.path.exists(f"{d}/species.tsv.gz")):
        return False, 0.0, 0  # rescue run whose world died before the partner was added back
    seq = species(d); last = None; pop = collections.Counter()
    for l in ropen(f"{d}/popdy001.dat"):
        t, s, c = l.split(","); t = int(t)
        if t != last:
            last = t; pop = collections.Counter()
        pop[seq.get(int(s), "")] += int(c)
    alive = "EXTINCT" not in open(f"{d}/log.txt").read()
    n = sum(pop.values())
    return alive, (sum(pop[x] for x in tracked) / n if n and alive else 0.0), sum(pop[x] for x in tracked) if alive else 0

def established(d, start):
    w = collections.defaultdict(collections.Counter)
    for l in ropen(f"{d}/copies.tsv"):
        r = l.rstrip("\n").split("\t")
        t = int(r[0])
        if t > start and r[2] == r[3] and len(r) > 6:
            w[(t - 1) // 10000][(int(r[1]), int(r[2]))] += int(r[6])
    best = 0.0
    for c in w.values():
        tot = sum(c.values())
        if not tot:
            continue
        mut = sum(v + c[(b, a)] for (a, b), v in c.items() if a < b and v >= 20 and c[(b, a)] >= 20)
        best = max(best, mut / tot)
    return best

def table(kind, tracked_of):
    print(f"\n== {kind}: condition, alive/runs, median final share of tracked sequences")
    for cond in sorted(os.listdir(f"{B2}/{kind}")):
        runs = sorted(d for d in glob.glob(f"{B2}/{kind}/{cond}/seed*") if os.path.exists(f"{d}/log.txt"))
        if not runs:
            continue
        res = [final(d, tracked_of(cond)) for d in runs]
        al = [r for r in res if r[0]]
        med = statistics.median(r[1] for r in al) if al else 0
        kept = sum(1 for r in al if r[2] > 0)
        print(f"  {cond:<14} alive {len(al):>2}/{len(res):<2}  tracked present {kept:>2}  median share {med:.2f}")

if os.path.isdir(f"{B2}/invade"):
    table("invade", lambda c: PAIRS[c.split("_")[0]])
if os.path.isdir(f"{B2}/reverse"):
    table("reverse", lambda c: (COP,))
if os.path.isdir(f"{B2}/bottleneck"):
    table("bottleneck", lambda c: PAIRS.get(c.split("_")[0], (COP,)))
if os.path.isdir(f"{B2}/rescue"):
    table("rescue", lambda c: PAIRS[c.split("_")[0]])
if os.path.isdir(f"{B2}/access"):
    print("\n== access: established partnership (best 10k-window share of placed exact copies in mutual pairs >= 0.10)")
    for grp in sorted(os.listdir(f"{B2}/access")):
        seed, t = grp.split("_t"); seed = seed.replace("seed", ""); t = int(t)
        for kind in ("replay", "intro_own", "intro_s1"):
            runs = sorted(d for d in glob.glob(f"{B2}/access/{grp}/{kind}*") if os.path.exists(f"{d}/log.txt"))
            if not runs:
                continue
            est = [established(d, t) for d in runs]
            line = f"  seed{seed} t={t//1000:>3}k {kind:<10} established {sum(e >= 0.10 for e in est)}/{len(runs)}  best shares {' '.join(f'{e:.2f}' for e in est)}"
            if kind != "replay":
                tr = OWN[seed] if kind == "intro_own" else PAIRS["s1"]
                fin = [final(d, tr) for d in runs]
                line += f"  | introduced present at end {sum(1 for f in fin if f[2] > 0)}/{len(runs)}"
            print(line)
