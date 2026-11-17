#!/usr/bin/env python3
"""compete.py setup OUT    -> builds the competition runs (jobs.txt)
   compete.py score OUT    -> scores them

Competition between a replicator without self-scan (A) and one with it (D), taken from the same
snapshot of a published run (the first 20k snapshot where a self-scanning replicator is among
the top replicators), with and without that snapshot's parasites, with mutation off and on.
Nothing in the chemistry is changed: only the starting mix differs.

  A, D        most abundant non-self-scanning / self-scanning self-copiers in the snapshot
  parasites   up to 5 commonest species in the snapshot's top 30 that cannot copy themselves but
              are copied exactly by A or D in the probe
  start       1000 A + 1000 D at random cells, + 1000 parasites (split by abundance) if added
  conditions  m0p0 (mutation off, no parasites), m0p1 (off, parasites), m1p0, m1p1 (MUTATE 0.0002)
  runs        10 seeds each, 30k steps, box-config rules (old rules + COMPAT_DELETION)
Score: share of D among molecules assigned to A or D. With mutation on, a species is assigned to
whichever of A or D it is nearer by edit distance, if within 20% of the length.
"""
import sys, os, re, random, collections, subprocess, tempfile, glob, statistics
H = os.path.expanduser("~")
RAW = f"{H}/published/bio/biopaperdata/raw"
REACT = f"{H}/smol-src/target/release/smol-react"
SMOL = f"{H}/smol-src/target/release/smol"
BASE = f"{H}/smol/restore/pub/box/seed1101/run.conf"
mode, OUT = sys.argv[1], sys.argv[2]
SNAPS = {1: 180000, 2: 200000, 3: 80000, 4: 100000, 5: 580000, 7: 80000}

def probe(seqs, flags, steps=5000):
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        for i, s in enumerate(seqs): f.write(f"R m{i} {s}\n")
        n = f.name
    o = subprocess.run([REACT, n, str(steps)] + flags, capture_output=True, text=True).stdout
    os.unlink(n); return o

def parse(out):
    res = {}; cur = None
    for l in out.splitlines():
        m = re.match(r"(m\d+)\s+initiates with (m\d+)\s+P\(bind\)=([\d.]+)", l)
        if m: cur = (int(m.group(1)[1:]), int(m.group(2)[1:])); res[cur] = set(); continue
        m = re.match(r"\s+-> (?:product|passive made) (m\d+)", l)
        if m and cur: res[cur].add(int(m.group(1)[1:]))
    return res

def expo(S):
    tr = probe([S], ["--steptrace", "--self"])
    nc = sum(1 for l in tr.splitlines() if l.startswith("T ") and l.split()[2] == "=")
    pl = sum(len(l.split()[2]) for l in tr.splitlines() if l.startswith("PROD"))
    return nc / pl if pl else 0

def snapshot(r, step):
    c = collections.Counter()
    for l in open(f"{RAW}/run{r}/out1_{step}.conf", errors="replace"):
        if l.startswith("AGENT"): c[l.split()[1]] += 1
        elif l.startswith("###"):
            p = l.split()
            if len(p) >= 3: c[p[2]] += 1
    return c

def lev(a, b):
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i] + [0] * len(b)
        for j, cb in enumerate(b, 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb))
        prev = cur
    return prev[-1]

if mode == "setup":
    os.makedirs(OUT, exist_ok=True)
    base = [l.rstrip("\n") for l in open(BASE) if l.split() and l.split()[0] not in
            ("RANDSEED", "NSTEPS", "OFFSPRINGLOG", "SNAPEVERY", "CKPTEVERY", "KEEPCKPT", "CKPTKEEPEVERY", "MUTATE", "AGENT", "GRIDPOS")]
    jobs = open(f"{OUT}/jobs.txt", "w"); meta = open(f"{OUT}/pairs.tsv", "w")
    for r, step in SNAPS.items():
        c = snapshot(r, step)
        top = [s for s, _ in c.most_common(30)]
        selfres = parse(probe(top, ["--self"]))
        selfcop = [s for i, s in enumerate(top) if i in selfres.get((i, i), set())]
        ex = {s: expo(s) for s in selfcop}
        A = next(s for s in selfcop if ex[s] < 1.5); D = next(s for s in selfcop if ex[s] >= 1.5)
        cand = [s for s in top if s not in selfcop][:20]
        res = parse(probe([A, D] + cand, []))
        par = [s for k, s in enumerate(cand, 2) if any(k in res.get(key, set()) for key in ((0, k), (k, 0), (1, k), (k, 1)))][:5]
        meta.write(f"run{r}\t{step}\t{A}\t{D}\t{','.join(par)}\t{c[A]}\t{c[D]}\n")
        for m in (0, 1):
            for p in (0, 1):
                if p and not par: continue
                for seed in range(1, 11):
                    d = f"{OUT}/run{r}/m{m}p{p}/s{seed}"; os.makedirs(d, exist_ok=True)
                    rng = random.Random(seed * 1000 + r)
                    cells = rng.sample([(x, y) for x in range(125) for y in range(100)], 3100)
                    mix = [(A, 1000), (D, 1000)]
                    if p:
                        tot = sum(c[s] for s in par)
                        mix += [(s, round(1000 * c[s] / tot)) for s in par]
                    lines = [f"RANDSEED {seed}", "NSTEPS 30000", f"MUTATE {0.0002 if m else 0}", "COMPAT_DELETION 1",
                             "SNAPEVERY 100000000", "CKPTEVERY 100000", "POPDYEVERY 500"] + base
                    k = 0
                    for s, n in mix:
                        if n <= 0: continue
                        lines.append(f"AGENT {s} {n} Q")
                        for _ in range(n):
                            x, y = cells[k]; k += 1
                            lines.append(f"GRIDPOS {x} {y}")
                    open(f"{d}/run.conf", "w").write("\n".join(lines) + "\n")
                    jobs.write(f"cd {d} && {SMOL} run run.conf > log.txt 2>&1\n")
    print(open(f"{OUT}/pairs.tsv").read())

else:
    pairs = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(f"{OUT}/pairs.tsv")}
    checkpoints = [2000, 5000, 10000, 20000, 30000]
    table = collections.defaultdict(list)
    for d in sorted(glob.glob(f"{OUT}/run*/m*p*/s*")):
        run, cond = d.split("/")[-3], d.split("/")[-2]
        A, D = pairs[run][2], pairs[run][3]
        seq = {}
        for l in open(f"{d}/species.tsv"):
            p = l.rstrip("\n").split("\t")
            if p[0].isdigit(): seq[int(p[0])] = p[-1]
        side = {}
        for sp, s in seq.items():
            if s == A: side[sp] = "A"
            elif s == D: side[sp] = "D"
            elif cond.startswith("m1"):
                da, dd = lev(s, A), lev(s, D)
                if min(da, dd) <= 0.2 * min(len(A), len(D)): side[sp] = "A" if da < dd else "D"
        cnt = collections.defaultdict(lambda: collections.Counter())
        for l in open(f"{d}/popdy001.dat"):
            st, sp, n = (int(v) for v in l.split(","))
            if st in checkpoints and sp in side: cnt[st][side[sp]] += n
        for st in checkpoints:
            a, dd = cnt[st]["A"], cnt[st]["D"]
            table[(run, cond, st)].append(dd / (a + dd) if a + dd else float("nan"))
    print("Share of self-scanners (D) among A+D; median over seeds [extinct or lost = nan excluded] (n)")
    conds = ["m0p0", "m0p1", "m1p0", "m1p1"]
    for st in checkpoints:
        print(f"== step {st}")
        for run in sorted(pairs, key=lambda x: int(x[3:])):
            row = []
            for c in conds:
                v = [x for x in table[(run, c, st)] if x == x]
                row.append(f"{c} {statistics.median(v):.2f} ({len(v)})" if v else f"{c}   -  ")
            print(f"  {run:5} " + "   ".join(row))
