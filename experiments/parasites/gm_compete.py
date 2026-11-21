#!/usr/bin/env python3
"""gm_compete.py setup OUT | probe OUT | score OUT

Can a change to the machine be selected for, and do gene-parasites make the difference?
(findings 34: viable machine variants exist and are linked to their genes, yet none spreads.)
Contest: the seed gene-machine pair (G0, M0) against a variant pair (Gv, Mv) whose machine has a
duplicated copy loop (Gv = G0's head + smol-codon encode(Mv)). Fresh 125x100 world: 500 M0 + 500 G0
+ 500 Mv + 500 Gv at random cells, with or without 1,000 gene-parasites (the 5 commonest class-P
species of gm run k at 1M, split by abundance; seed k uses run 30k's pool), mutation off or on,
10 seeds, 30k steps, the gm runs' rules. Score: share of Gv among Gv + G0 at 30k (with mutation on,
genes within 10 edits are assigned to the nearer of the two)."""
import sys, os, random, glob, statistics, subprocess, collections, re, tempfile
H = os.path.expanduser("~"); B = f"{H}/smol-src/target/release"
mode, OUT = sys.argv[1], os.path.abspath(sys.argv[2])
M0 = "WGEWLHHRJXQQ^B?WGEW}$=?>%^F$LUBO>B^C$=?>%}"
G0 = "JTRJYUUEWKQQYHOBAOAFAIAOACADADCVALBIAJAJBGAXBQAOAFAIAOFWCZAABQAPBZBGCPCZACBEAXABAPAXBGAMCZAABQAPBZFW"
VARS = {"loop1": "WGEWLHHRJXQQ^B?WGEW}$=?>$=?>%^F$LUBO>B^C$=?>%}",
        "loop2": "WGEWLHHRJXQQ^B?WGEW}$=?>%^F$LUBO>B^C$=?>$=?>%}",
        "half": "WGEWLHHRJXQQ^B?WGEW}$=?$=?>%^F$LUBO>B^C$=?>%}"}
def gene(m): return G0[:16] + subprocess.run([f"{B}/smol-codon", "encode", m], capture_output=True, text=True).stdout.strip()
BASE = [l for l in open(f"{H}/smol/runs/gm/s301/run.conf").read().splitlines()
        if l.split() and l.split()[0] not in ("RANDSEED", "NSTEPS", "MUTATE", "AGENT", "GRIDPOS", "CKPTEVERY", "KEEPCKPT", "CKPTKEEPEVERY")]
def pool(k):
    w = f"{H}/smol/runs/ts_gm/s{300 + k}/w1000k"
    seq = {}
    for l in open(f"{w}/species.tsv"):
        p = l.rstrip("\n").split("\t")
        if p[0].isdigit() and len(p) >= 5: seq.setdefault(int(p[0]), p[-1])
    cl = [l.split() for l in open(f"{w}/classes.tsv")]
    ps = [(int(n), seq[int(s)]) for s, c, n in cl if c == "P" and int(s) in seq][:5]
    return [s for _, s in ps], [n for n, _ in ps]
def lev(a, b):
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i] + [0] * len(b)
        for j, cb in enumerate(b, 1): cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb))
        prev = cur
    return prev[-1]

if mode == "setup":
    os.makedirs(OUT, exist_ok=True); jobs = open(f"{OUT}/jobs.txt", "w")
    for v, Mv in VARS.items():
        Gv = gene(Mv)
        for mut in (0, 1):
            for p in (0, 1):
                for seed in range(1, 11):
                    d = f"{OUT}/{v}/m{mut}p{p}/s{seed}"; os.makedirs(d, exist_ok=True)
                    mix = [(M0, 500), (G0, 500), (Mv, 500), (Gv, 500)]
                    if p:
                        ps, ns = pool(seed); t = sum(ns)
                        mix += [(s, round(1000 * n / t)) for s, n in zip(ps, ns)]
                    rng = random.Random(seed * 1000 + 7)
                    cells = rng.sample([(x, y) for x in range(125) for y in range(100)], sum(n for _, n in mix))
                    lines = [f"RANDSEED {seed}", "NSTEPS 30000", f"MUTATE {0.0002 if mut else 0}"] + [l for l in BASE if not l.startswith(("POPDYEVERY", "COPYEVERY", "SNAPEVERY"))] + \
                            ["POPDYEVERY 500", "SNAPEVERY 100000000", "CKPTEVERY 1000000"]
                    k = 0
                    for s, n in mix:
                        if n <= 0: continue
                        lines.append(f"AGENT {s} {n} Q")
                        for _ in range(n):
                            x, y = cells[k]; k += 1; lines.append(f"GRIDPOS {x} {y}")
                    open(f"{d}/run.conf", "w").write("\n".join(lines) + "\n")
                    jobs.write(f"cd {d} && {B}/smol run run.conf > log.txt 2>&1; rm -f copies.tsv grid.snap\n")
    jobs.close(); print(open(f"{OUT}/jobs.txt").read().count("\n"), "jobs")

elif mode == "probe":
    print("machine vs each run's 5 gene-parasites (abundance-weighted): share of bound time copied, mean steps bound; self-loop steps")
    for name, M in [("seed", M0)] + list(VARS.items()):
        G = gene(M) if name != "seed" else G0
        cop = []; stp = []
        for k in range(1, 11):
            ps, ns = pool(k)
            with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
                f.write(f"R M {M}\n" + "".join(f"R p{i} {s}\n" for i, s in enumerate(ps))); fn = f.name
            out = subprocess.run([f"{B}/smol-react", fn, "5000", "--strand", "--translate", "1", "--stop", "--hub"], capture_output=True, text=True).stdout
            os.unlink(fn)
            tc = tt = st = 0.0
            for i, n in enumerate(ns):
                m = re.search(rf"^M\s+initiates with p{i}\s+P\(bind\)=([\d.]+).*?steps=(\d+)(.*?)(?=^\S|\Z)", out, re.M | re.S)
                if not m: continue
                pb, s = float(m.group(1)), int(m.group(2)); w = n * pb * s
                tt += w; st += n * pb * s
                if re.search(rf"-> (?:product|passive made) p{i}\b", m.group(3)): tc += w
            if tt: cop.append(tc / tt); stp.append(st / sum(ns))
        print(f"  {name:6} copied {statistics.median(cop):.2f}  bound-steps x P(bind) {statistics.median(stp):6.0f}")

else:
    print("share of the variant's gene among both genes at 30k: median over 10 seeds (seeds where both genes are gone excluded)")
    for v, Mv in VARS.items():
        Gv = gene(Mv); row = []
        for c in ("m0p0", "m0p1", "m1p0", "m1p1"):
            vals = []
            for d in sorted(glob.glob(f"{OUT}/{v}/{c}/s*")):
                seq = {}
                for l in open(f"{d}/species.tsv"):
                    p = l.rstrip("\n").split("\t")
                    if p[0].isdigit() and len(p) >= 5: seq[int(p[0])] = p[-1]
                alive = collections.Counter()
                for l in open(f"{d}/popdy001.dat"):
                    st, sp, n = (int(x) for x in l.split(","))
                    if st == 30000: alive[sp] += n
                side = {}
                for sp, s in ((sp, seq[sp]) for sp in alive if sp in seq):
                    if s == Gv: side[sp] = "v"
                    elif s == G0: side[sp] = "0"
                    elif c.startswith("m1") and sum(ch in "$>^?=%}" for ch in s) <= 2 and abs(len(s) - len(G0)) < 16:
                        a, b = lev(s, G0), lev(s, Gv)
                        if min(a, b) <= 10 and a != b: side[sp] = "v" if b < a else "0"
                cnt = collections.Counter()
                for sp, n in alive.items():
                    if sp in side: cnt[side[sp]] += n
                if cnt["v"] + cnt["0"]: vals.append(cnt["v"] / (cnt["v"] + cnt["0"]))
            row.append(f"{c} {statistics.median(vals):.2f} ({len(vals)})" if vals else f"{c}   -  ")
        print(f"  {v:6} {Mv}\n         " + "   ".join(row))
