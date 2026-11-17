#!/usr/bin/env python3
"""knockout.py setup OUT  /  knockout.py score OUT

Pools where the fastest-parasite rule (predict.py E3) gets the winner wrong: is it because some
pool members are themselves copiers (molecules that cannot copy themselves but copy other
parasites), giving a parasite extra hosts besides A or D? Each such pool is rerun with its copier
members removed (the rest rescaled to 1000), same seeds and set-up as before. Controls: four pools
the rule got right that also contain copier members, rerun the same way.
A pool member is a copier if, in a pairwise probe (smol-react --compat) with A, D and the rest of
the pool, it makes an exact copy of any other molecule in the set.
"""
import sys, os, re, glob, random, statistics, subprocess, tempfile
H = os.path.expanduser("~")
mode, OUT = sys.argv[1], os.path.abspath(sys.argv[2])
R = f"{H}/smol-src/target/release/smol-react"; SMOL = f"{H}/smol-src/target/release/smol"
PAIRS = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(f"{H}/smol/runs/compete/pairs.tsv")}
BASECONF = f"{H}/smol/restore/pub/box/seed1101/run.conf"
BASE = [l.rstrip("\n") for l in open(BASECONF) if l.split() and l.split()[0] not in
        ("RANDSEED", "NSTEPS", "OFFSPRINGLOG", "SNAPEVERY", "CKPTEVERY", "KEEPCKPT", "CKPTKEEPEVERY", "MUTATE", "AGENT", "GRIDPOS", "POPDYEVERY")]

def pools():
    out = {}
    for l in open(f"{H}/smol/runs/ratio/chosen.tsv"):
        q = l.rstrip("\n").split("\t")
        out[(q[0], f"src{q[1]}")] = (q[6].split(","), [int(x) for x in q[7].split(",")], f"{H}/smol/runs/ratio/{q[0]}/src{q[1]}")
    for l in open(f"{H}/smol/runs/evolve/pools.tsv"):
        q = l.rstrip("\n").split("\t")
        if int(q[5]) < 200 or not q[8]: continue
        nm = f"on{q[1]}{q[2]}_{int(q[3]) // 1000}k"
        out[(q[0], nm)] = (q[8].split(","), [int(x) for x in q[9].split(",")], f"{H}/smol/runs/evolve/{q[0]}/comp{q[1]}{q[2]}_{int(q[3]) // 1000}k")
    return out

def copiers(A, D, par):
    seqs = [A, D] + par
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        for i, s in enumerate(seqs): f.write(f"R m{i} {s}\n")
        n = f.name
    out = subprocess.run([R, n, "5000", "--compat"], capture_output=True, text=True).stdout; os.unlink(n)
    cop = set(); cur = None
    for l in out.splitlines():
        m = re.match(r"m\d+\s+initiates with m\d+\s+P\(bind\)=[\d.]+\s+copier=m(\d+)\s+copied=m(\d+)", l)
        if m: cur = (int(m.group(1)), int(m.group(2))); continue
        m = re.match(r"\s+-> (?:product|passive made) m(\d+)", l)
        if m and cur and int(m.group(1)) == cur[1] and cur[0] >= 2 and cur[0] != cur[1]: cop.add(cur[0] - 2)
    return cop

def share(d, A, D):
    seq = {}
    for l in open(f"{d}/species.tsv"):
        p = l.rstrip("\n").split("\t")
        if p[0].isdigit(): seq[int(p[0])] = p[-1]
    a = b = 0
    for l in open(f"{d}/popdy001.dat"):
        s, sp, n = (int(v) for v in l.split(","))
        if s == 30000 and sp in seq:
            if seq[sp] == A: a += n
            elif seq[sp] == D: b += n
    return b / (a + b) if a + b else None

def med(ds, A, D):
    v = [x for x in (share(d, A, D) for d in ds) if x is not None]
    return statistics.median(v) if v else float("nan")

P = pools()
pred = {}
for l in open(f"{H}/smol/runs/predict2.txt"):
    q = l.split()
    if len(q) > 12 and q[0] in ("snap", "evo"):
        pred[(q[1], q[2])] = (float(q[7]), float(q[8]), q[-1])   # E3(A), E3(D), ok/MISS

if mode == "setup":
    os.makedirs(OUT, exist_ok=True)
    jobs = open(f"{OUT}/jobs.txt", "w"); info = open(f"{OUT}/pools.tsv", "w")
    misses = [k for k, v in pred.items() if v[2] == "MISS"]
    ctl = []
    for k, v in sorted(pred.items()):
        if v[2] == "ok" and len(ctl) < 4 and k in P:
            par = P[k][0]; A, D = PAIRS[k[0]][2], PAIRS[k[0]][3]
            if copiers(A, D, par) and k[0] not in [c[0] for c in ctl]: ctl.append(k)
    for kind, keys in (("miss", misses), ("ctl", ctl)):
        for run, nm in keys:
            par, cnt, _ = P[(run, nm)]; A, D = PAIRS[run][2], PAIRS[run][3]
            cop = copiers(A, D, par)
            keep = [(s, n) for i, (s, n) in enumerate(zip(par, cnt)) if i not in cop]
            info.write("\t".join([kind, run, nm, ",".join(par[i] for i in sorted(cop)), ",".join(s for s, _ in keep)]) + "\n")
            print(kind, run, nm, "copier members removed:", len(cop), "of", len(par), flush=True)
            if not cop or not keep: continue
            t = sum(n for _, n in keep)
            for cs in range(1, 11):
                d = f"{OUT}/{run}/{nm}/s{cs}"; os.makedirs(d, exist_ok=True)
                rng = random.Random(cs * 1000 + int(run[3:]))
                mix = [(A, 1000), (D, 1000)] + [(s, round(1000 * n / t)) for s, n in keep]
                cells = rng.sample([(x, y) for x in range(125) for y in range(100)], sum(n for _, n in mix))
                lines = [f"RANDSEED {cs}", "NSTEPS 30000", "MUTATE 0", "COMPAT_DELETION 1",
                         "SNAPEVERY 100000000", "CKPTEVERY 1000000", "POPDYEVERY 500"] + BASE
                k = 0
                for s, n in mix:
                    if n <= 0: continue
                    lines.append(f"AGENT {s} {n} Q")
                    for _ in range(n):
                        x, y = cells[k]; k += 1; lines.append(f"GRIDPOS {x} {y}")
                open(f"{d}/run.conf", "w").write("\n".join(lines) + "\n")
                jobs.write(f"cd {d} && {SMOL} run run.conf > log.txt 2>&1\n")
elif mode == "score":
    base = {run: med(glob.glob(f"{H}/smol/runs/compete/{run}/m0p0/s*"), p[2], p[3]) for run, p in PAIRS.items()}
    print("kind  run   pool         copiers  E3(A) E3(D) rule says   D share: full pool -> copiers removed   base   rule right after?")
    for l in open(f"{OUT}/pools.tsv"):
        kind, run, nm, cop, keep = l.rstrip("\n").split("\t")
        A, D = PAIRS[run][2], PAIRS[run][3]
        ea, ed, _ = pred[(run, nm)]
        says = "D wins" if ed < ea else "D loses"
        full = med(glob.glob(f"{P[(run, nm)][2]}/s*"), A, D)
        ko_d = glob.glob(f"{OUT}/{run}/{nm}/s*")
        ko = med(ko_d, A, D) if ko_d else float("nan")
        right = "-" if ko != ko else ("yes" if (ko > base[run]) == (ed < ea) else "no")
        print(f"{kind:5} {run:5} {nm:12} {len([c for c in cop.split(',') if c]):>3}     {ea:5.2f} {ed:5.2f} {says:9}  {full:5.2f} -> {ko:5.2f}                        {base[run]:.2f}   {right}")
