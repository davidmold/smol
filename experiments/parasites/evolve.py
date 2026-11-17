#!/usr/bin/env python3
"""evolve.py evo OUT      -> evolution worlds: 2000 copies of A alone or D alone, MUTATE 0.0002,
                            150k steps, 3 seeds, for each competition run (jobs_evo.txt)
   evolve.py comp OUT     -> harvests a parasite pool from each evolution world at 50k and 150k
                            and sets up A-vs-D competitions against it (jobs_comp.txt, pools.tsv)
   evolve.py score OUT    -> D share at 30k by run, host and harvest time

Causal test of "parasites adapt to the scanner": parasites evolved in a D-only world should beat
D in competition; parasites evolved in an A-only world should not. Same A, D, pool selection,
competition set-up and baseline as ratio.py / compete.py.
"""
import sys, os, collections, random, glob, statistics, math, types
H = os.path.expanduser("~")
mode, OUT = sys.argv[1], sys.argv[2]
C = types.SimpleNamespace()
_src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "compete.py")).read()
_src = _src.split('if mode == "setup":')[0].replace("mode, OUT = sys.argv[1], sys.argv[2]", "")
exec(_src, C.__dict__)
PAIRS = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(f"{H}/smol/runs/compete/pairs.tsv")}
BASE = [l.rstrip("\n") for l in open(C.BASE) if l.split() and l.split()[0] not in
        ("RANDSEED", "NSTEPS", "OFFSPRINGLOG", "SNAPEVERY", "CKPTEVERY", "KEEPCKPT", "CKPTKEEPEVERY", "MUTATE", "AGENT", "GRIDPOS", "POPDYEVERY")]
HARVEST = (50000, 150000)
SEEDS = (1, 2, 3)

def conf(seed, nsteps, mutate, mix, rngkey):
    rng = random.Random(rngkey)
    cells = rng.sample([(x, y) for x in range(125) for y in range(100)], sum(n for _, n in mix if n > 0))
    lines = [f"RANDSEED {seed}", f"NSTEPS {nsteps}", f"MUTATE {mutate}", "COMPAT_DELETION 1",
             "SNAPEVERY 100000000", "CKPTEVERY 1000000", "POPDYEVERY 500"] + BASE
    k = 0
    for s, n in mix:
        if n <= 0: continue
        lines.append(f"AGENT {s} {n} Q")
        for _ in range(n):
            x, y = cells[k]; k += 1
            lines.append(f"GRIDPOS {x} {y}")
    return "\n".join(lines) + "\n"

def counts_at(d, step):
    seq = {}
    for l in open(f"{d}/species.tsv"):
        p = l.rstrip("\n").split("\t")
        if p[0].isdigit(): seq[int(p[0])] = p[-1]
    c = collections.Counter(); last = collections.Counter(); laststep = 0
    for l in open(f"{d}/popdy001.dat"):
        s, sp, n = (int(v) for v in l.split(","))
        if s > step: break
        if s != laststep: last, laststep = collections.Counter(), s
        last[seq.get(sp, "")] += n
    return laststep, last

def pool(c, A, D):
    top = [s for s, _ in c.most_common(30) if s]
    selfres = C.parse(C.probe(top, ["--self"]))
    selfcop = {s for i, s in enumerate(top) if i in selfres.get((i, i), set())}
    cand = [s for s in top if s not in selfcop and s not in (A, D)][:20]
    if not cand: return []
    res = C.parse(C.probe([A, D] + cand, []))
    return [s for k, s in enumerate(cand, 2)
            if any(k in res.get(key, set()) for key in ((0, k), (k, 0), (1, k), (k, 1)))][:5]

os.makedirs(OUT, exist_ok=True)
if mode == "evo":
    jobs = open(f"{OUT}/jobs_evo.txt", "w")
    for run, p in PAIRS.items():
        A, D = p[2], p[3]
        for host, S in (("A", A), ("D", D)):
            for seed in SEEDS:
                d = f"{OUT}/{run}/evo{host}/e{seed}"; os.makedirs(d, exist_ok=True)
                open(f"{d}/run.conf", "w").write(conf(100 + seed, max(HARVEST), 0.0002, [(S, 2000)], seed * 7919 + int(run[3:])))
                jobs.write(f"cd {d} && {C.SMOL} run run.conf > log.txt 2>&1\n")

elif mode == "comp":
    jobs = open(f"{OUT}/jobs_comp.txt", "w"); pf = open(f"{OUT}/pools.tsv", "w")
    for run, p in PAIRS.items():
        A, D = p[2], p[3]
        for host in "AD":
            for seed in SEEDS:
                d = f"{OUT}/{run}/evo{host}/e{seed}"
                for hv in HARVEST:
                    st, c = counts_at(d, hv)
                    tot = sum(c.values())
                    par = pool(c, A, D) if tot else []
                    cnt = [c[s] for s in par]
                    pf.write("\t".join(map(str, [run, host, seed, hv, st, tot, c[A], c[D], ",".join(par), ",".join(map(str, cnt))])) + "\n"); pf.flush()
                    print(run, host, seed, hv, "step", st, "pop", tot, "A", c[A], "D", c[D], "parasites", len(par), flush=True)
                    if not par: continue
                    t = sum(cnt)
                    for cs in range(1, 11):
                        cd = f"{OUT}/{run}/comp{host}{seed}_{hv // 1000}k/s{cs}"; os.makedirs(cd, exist_ok=True)
                        mix = [(A, 1000), (D, 1000)] + [(s, round(1000 * n / t)) for s, n in zip(par, cnt)]
                        open(f"{cd}/run.conf", "w").write(conf(cs, 30000, 0, mix, cs * 1000 + int(run[3:])))
                        jobs.write(f"cd {cd} && {C.SMOL} run run.conf > log.txt 2>&1\n")

else:
    def share(d, A, D, st=30000):
        seq = {}
        for l in open(f"{d}/species.tsv"):
            p = l.rstrip("\n").split("\t")
            if p[0].isdigit(): seq[int(p[0])] = p[-1]
        a = dd = 0
        for l in open(f"{d}/popdy001.dat"):
            s, sp, n = (int(v) for v in l.split(","))
            if s == st and sp in seq:
                if seq[sp] == A: a += n
                elif seq[sp] == D: dd += n
        return a, dd
    med = lambda v: statistics.median(v) if v else float("nan")
    print("D share of A+D at 30k (median of 10 seeds); A/D = seeds where that replicator is gone; base = no-parasite baseline")
    summary = collections.defaultdict(list)
    poppool = {}
    for l in open(f"{OUT}/pools.tsv"):
        q = l.rstrip("\n").split("\t"); poppool[(q[0], q[1], int(q[2]), int(q[3]))] = int(q[5])
    for run in sorted(PAIRS, key=lambda x: int(x[3:])):
        A, D = PAIRS[run][2], PAIRS[run][3]
        bl = med([b / (a + b) for a, b in (share(d, A, D) for d in glob.glob(f"{H}/smol/runs/compete/{run}/m0p0/s*")) if a + b])
        print(f"== {run}  baseline {bl:.2f}")
        for hv in HARVEST:
            row = []
            for host in "AD":
                for seed in SEEDS:
                    ds = glob.glob(f"{OUT}/{run}/comp{host}{seed}_{hv // 1000}k/s*")
                    if poppool.get((run, host, seed, hv), 0) < 200: row.append(f"{host}{seed} dead "); continue
                    if not ds: row.append(f"{host}{seed}  -  "); continue
                    v = [share(d, A, D) for d in ds]
                    sh = [b / (a + b) for a, b in v if a + b]
                    ga, gd = sum(1 for a, b in v if a == 0), sum(1 for a, b in v if b == 0)
                    m = med(sh); summary[(host, hv)].append((m, bl))
                    row.append(f"{host}{seed} {m:.2f} ({ga}/{gd})")
            print(f"  {hv // 1000:>3}k  " + "   ".join(row))
    print("\nOver runs x evolution seeds:")
    for (host, hv), v in sorted(summary.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        v2 = [(m, b) for m, b in v if m == m]
        print(f"  parasites evolved on {host}, harvested {hv // 1000:>3}k: median D share {med([m for m, b in v2]):.2f}, "
              f"median minus baseline {med([m - b for m, b in v2]):+.2f}, D below baseline in "
              f"{sum(1 for m, b in v2 if m < b)}/{len(v2)}, D share < 0.1 in {sum(1 for m, b in v2 if m < 0.1)}/{len(v2)}")
