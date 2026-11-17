#!/usr/bin/env python3
"""ratio.py survey OUT  -> for each run, D:A ratio and parasite pool in snapshots around the
                          competition snapshot (survey.tsv)
   ratio.py setup OUT   -> picks up to 5 snapshots per run spanning the D:A ratio, writes jobs.txt
   ratio.py score OUT   -> D share at 10k/20k/30k against the source snapshot's D:A ratio

Do parasites taken from where D was dominant beat D more often? Same A and D as compete.py
(compete/pairs.tsv); only the parasite pool changes, drawn from a different snapshot of the same
published run. Parasites chosen exactly as in compete.py: up to 5 commonest species in the
snapshot's top 30 that cannot copy themselves but are copied exactly by A or D in the probe.
Ratio: molecules in the snapshot within 20% edit distance of D vs of A (nearest wins).
Competition: 1000 A + 1000 D + 1000 parasites (split by abundance), mutation off, 10 seeds, 30k
steps, box rules + COMPAT_DELETION. Baseline without parasites = compete m0p0 (same A, D, seeds).
"""
import sys, os, collections, random, glob, statistics, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
H = os.path.expanduser("~")
mode, OUT = sys.argv[1], sys.argv[2]
import types
C = types.SimpleNamespace()
_src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "compete.py")).read()
_src = _src.split('if mode == "setup":')[0].replace("mode, OUT = sys.argv[1], sys.argv[2]", "")
exec(_src, C.__dict__)

PAIRS = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(f"{H}/smol/runs/compete/pairs.tsv")}

def family(c, A, D):
    a = d = 0
    lim = 0.2 * min(len(A), len(D))
    for s, n in c.items():
        if abs(len(s) - len(A)) > lim and abs(len(s) - len(D)) > lim: continue
        da, dd = C.lev(s, A), C.lev(s, D)
        if min(da, dd) <= lim:
            if dd < da: d += n
            else: a += n
    return a, d

def pool(c, A, D):
    top = [s for s, _ in c.most_common(30)]
    selfres = C.parse(C.probe(top, ["--self"]))
    selfcop = {s for i, s in enumerate(top) if i in selfres.get((i, i), set())}
    cand = [s for s in top if s not in selfcop and s not in (A, D)][:20]
    if not cand: return []
    res = C.parse(C.probe([A, D] + cand, []))
    return [s for k, s in enumerate(cand, 2)
            if any(k in res.get(key, set()) for key in ((0, k), (k, 0), (1, k), (k, 1)))][:5]

os.makedirs(OUT, exist_ok=True)
if mode == "survey":
    f = open(f"{OUT}/survey.tsv", "w")
    for run, p in PAIRS.items():
        r, snap, A, D = int(run[3:]), int(p[1]), p[2], p[3]
        for step in range(max(20000, snap - 60000), snap + 220001, 20000):
            try: c = C.snapshot(r, step)
            except FileNotFoundError: continue
            a, d = family(c, A, D)
            par = pool(c, A, D)
            tot = sum(c[s] for s in par)
            line = [run, step, c[A], c[D], a, d, ",".join(par), ",".join(str(c[s]) for s in par), tot]
            f.write("\t".join(map(str, line)) + "\n"); f.flush()
            print(*line[:6], len(par), tot, flush=True)

elif mode == "setup":
    rows = collections.defaultdict(list)
    for l in open(f"{OUT}/survey.tsv"):
        q = l.rstrip("\n").split("\t")
        if q[6] and int(q[4]) + int(q[5]) > 0: rows[q[0]].append(q)
    base = [l.rstrip("\n") for l in open(C.BASE) if l.split() and l.split()[0] not in
            ("RANDSEED", "NSTEPS", "OFFSPRINGLOG", "SNAPEVERY", "CKPTEVERY", "KEEPCKPT", "CKPTKEEPEVERY", "MUTATE", "AGENT", "GRIDPOS")]
    jobs = open(f"{OUT}/jobs.txt", "w"); chosen = open(f"{OUT}/chosen.tsv", "w")
    for run, qs in rows.items():
        lr = lambda q: math.log((int(q[5]) + 1) / (int(q[4]) + 1))
        qs.sort(key=lr)
        orig = [q for q in qs if int(q[1]) == int(PAIRS[run][1])]
        idx = sorted({0, len(qs) - 1, len(qs) // 4, len(qs) // 2, 3 * len(qs) // 4})
        pick = {q[1]: q for q in [qs[i] for i in idx] + orig}
        A, D = PAIRS[run][2], PAIRS[run][3]
        for step, q in sorted(pick.items(), key=lambda kv: int(kv[0])):
            chosen.write("\t".join(q) + "\n")
            par = q[6].split(","); cnt = [int(x) for x in q[7].split(",")]; tot = sum(cnt)
            for seed in range(1, 11):
                d = f"{OUT}/{run}/src{step}/s{seed}"; os.makedirs(d, exist_ok=True)
                rng = random.Random(seed * 1000 + int(run[3:]))
                cells = rng.sample([(x, y) for x in range(125) for y in range(100)], 3100)
                mix = [(A, 1000), (D, 1000)] + [(s, round(1000 * n / tot)) for s, n in zip(par, cnt)]
                lines = [f"RANDSEED {seed}", "NSTEPS 30000", "MUTATE 0", "COMPAT_DELETION 1",
                         "SNAPEVERY 100000000", "CKPTEVERY 100000", "POPDYEVERY 500"] + base
                k = 0
                for s, n in mix:
                    if n <= 0: continue
                    lines.append(f"AGENT {s} {n} Q")
                    for _ in range(n):
                        x, y = cells[k]; k += 1
                        lines.append(f"GRIDPOS {x} {y}")
                open(f"{d}/run.conf", "w").write("\n".join(lines) + "\n")
                jobs.write(f"cd {d} && {C.SMOL} run run.conf > log.txt 2>&1\n")
    print(open(f"{OUT}/chosen.tsv").read())

else:
    def share(d, A, D, st):
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
        return dd / (a + dd) if a + dd else float("nan")
    med = lambda v: statistics.median([x for x in v if x == x]) if any(x == x for x in v) else float("nan")
    print("run    source   A-fam  D-fam  log(D/A)  parasites  | D share at 10k 20k 30k (median of 10)   D extinct by 30k   | no-parasite baseline 30k")
    allpts = []
    for l in open(f"{OUT}/chosen.tsv"):
        q = l.rstrip("\n").split("\t"); run, step = q[0], q[1]
        A, D = PAIRS[run][2], PAIRS[run][3]
        ds = sorted(glob.glob(f"{OUT}/{run}/src{step}/s*"))
        sh = {st: [share(d, A, D, st) for d in ds] for st in (10000, 20000, 30000)}
        ext = sum(1 for x in sh[30000] if x == 0.0)
        bl = med([share(d, A, D, 30000) for d in glob.glob(f"{H}/smol/runs/compete/{run}/m0p0/s*")])
        lr = math.log((int(q[5]) + 1) / (int(q[4]) + 1))
        tag = " *" if step == PAIRS[run][1] else ""
        print(f"{run:6} {step:>7}{tag:2} {q[4]:>6} {q[5]:>6} {lr:8.2f}  {len(q[6].split(',')):>3}        |  "
              f"{med(sh[10000]):.2f} {med(sh[20000]):.2f} {med(sh[30000]):.2f}            {ext:>2}/{len(ds)}               |  {bl:.2f}")
        allpts.append((run, lr, med(sh[30000]) - bl))
    print("\n* = snapshot used in the original competition")
    # within-run rank correlation of log(D/A) with D's 30k share relative to baseline
    for run in sorted({p[0] for p in allpts}, key=lambda x: int(x[3:])):
        pts = [(x, y) for r, x, y in allpts if r == run and y == y]
        if len(pts) < 3: continue
        rk = lambda v: [sorted(v).index(x) for x in v]
        xs, ys = rk([p[0] for p in pts]), rk([p[1] for p in pts])
        n = len(pts); rho = 1 - 6 * sum((a - b) ** 2 for a, b in zip(xs, ys)) / (n * (n * n - 1))
        print(f"{run}: Spearman rho(log D/A, D share - baseline) = {rho:+.2f} over {n} pools")
