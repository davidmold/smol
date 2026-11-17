#!/usr/bin/env python3
"""dip.py RUN_DIR... : are parasites like defective interfering particles (regenerated from the
replicator each time) or like an evolving parasite lineage?

For each run (needs the full species.tsv with parents, copies.tsv, checkpoints every 100k):
  classes   over the whole run: R made exact copies; P was copied exactly >= 20 times and never
            made a copy (a successful parasite)
  origin    each species' template parent is whichever of its two reactant parents it is nearer
            to (edit distance); a parasite's lineage is followed back through P parents until a
            non-P ancestor. depth 1 = made directly from a replicator (regenerated); depth >= 2 =
            descended from another parasite
  exploit   in the probe, against the top replicator of the next 100k checkpoint: P(bind) x
            (the replicator copies it). Compared between depth 1 and depth >= 2 parasites
  age       at each 100k checkpoint: median age of the parasite molecules' lineages (from the
            founding parasite's birth) and of the replicator molecules' species
"""
import sys, os, re, gzip, glob, collections, statistics, subprocess, tempfile, random

REACT = os.path.expanduser("~/smol-src/target/release/smol-react")

def lev(a, b):
    if a == b: return 0
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i] + [0] * len(b)
        for j, cb in enumerate(b, 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb))
        prev = cur
    return prev[-1]

def probe(seqs):
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        for i, s in enumerate(seqs): f.write(f"R m{i} {s}\n")
        n = f.name
    out = subprocess.run([REACT, n, "5000", "--hub"], capture_output=True, text=True).stdout
    os.unlink(n)
    res = {}; cur = None
    for l in out.splitlines():
        m = re.match(r"(m\d+)\s+initiates with (m\d+)\s+P\(bind\)=([\d.]+)", l)
        if m: cur = (int(m.group(1)[1:]), int(m.group(2)[1:])); res[cur] = [float(m.group(3)), set()]; continue
        m = re.match(r"\s+-> (?:product|passive made) (m\d+)", l)
        if m and cur: res[cur][1].add(int(m.group(1)[1:]))
    return res

def snapshot(path):
    cnt = collections.Counter()
    with gzip.open(path, "rt") as f:
        for l in f:
            if l.startswith("MOLS"): break
        for l in f: cnt[int(l.split()[3])] += 1
    return cnt

summary = []
for d in sys.argv[1:]:
    name = "/".join(d.rstrip("/").split("/")[-2:])
    seq, born, par = {}, {}, {}
    for l in open(f"{d}/species.tsv"):
        p = l.rstrip("\n").split("\t")
        if not p[0].isdigit(): continue
        i = int(p[0]); seq[i] = p[4]; born[i] = int(p[1]); par[i] = (int(p[2]), int(p[3]))
    copier, copied = set(), collections.Counter()
    for l in open(f"{d}/copies.tsv"):
        p = l.split("\t")
        if len(p) < 5 or not p[0].isdigit(): continue
        a, pa, c, n = int(p[1]), int(p[2]), int(p[3]), int(p[4])
        if c == pa or c == a:
            copier.add(a); copied[c] += n
    P = {s for s, n in copied.items() if n >= 20 and s not in copier}
    R = copier
    tmpl = {}
    def template(s):
        if s not in tmpl:
            a, b = par.get(s, (-1, -1))
            cands = [x for x in (a, b) if x in seq]
            tmpl[s] = min(cands, key=lambda x: lev(seq[s], seq[x])) if cands else -1
        return tmpl[s]
    depth, origin, founder = {}, {}, {}
    for s in P:
        k, cur = 1, s
        while True:
            t = template(cur)
            if t in P and k < 500:
                k += 1; cur = t; continue
            origin[s] = "R" if t in R else ("seed" if t == -1 else "other")
            founder[s] = cur; break
        depth[s] = k
    dcount = collections.Counter(min(v, 5) for v in depth.values())
    wtot = sum(copied[s] for s in P) or 1
    wd1 = sum(copied[s] for s in P if depth[s] == 1) / wtot
    ocount = collections.Counter(origin.values())
    clades = collections.Counter(founder.values())
    # exploitation vs depth, per checkpoint
    expl = {1: [], 2: []}
    cks = sorted(glob.glob(f"{d}/checkpoint_*00000.txt.gz"), key=lambda p: int(re.search(r"_(\d+)\.", p).group(1)))
    rng = random.Random(1)
    ages_p, ages_r, sages_p = [], [], []
    for ck in cks:
        st = int(re.search(r"_(\d+)\.", ck).group(1))
        cnt = snapshot(ck)
        reps = [s for s, _ in cnt.most_common() if s in R]
        if not reps: continue
        S = reps[0]
        pm = [s for s in cnt if s in P]
        if pm:
            ages_p.append(statistics.median([st - born[founder[s]] for s in pm for _ in range(min(cnt[s], 50))]))
            sages_p.append(statistics.median([st - born[s] for s in pm for _ in range(min(cnt[s], 50))]))
        rm = [s for s in cnt if s in R]
        ages_r.append(statistics.median([st - born[s] for s in rm for _ in range(min(cnt[s], 50))]))
        recent = [s for s in P if st - 100000 < born[s] <= st]
        rng.shuffle(recent); recent = recent[:150]
        if not recent: continue
        res = probe([seq[S]] + [seq[s] for s in recent])
        for j, s in enumerate(recent, 1):
            pb, prods = res.get((0, j), [0.0, set()])
            expl[1 if depth[s] == 1 else 2].append(pb * (j in prods))
    m = lambda v: statistics.mean(v) if v else float("nan")
    line = (f"{name:14} parasites {len(P):5d} (clades {len(clades):4d}; largest clade {max(clades.values()) if clades else 0:3d} species)  "
            f"depth1 {dcount[1]/max(len(P),1):.2f} depth2 {dcount[2]/max(len(P),1):.2f} depth3+ {sum(dcount[k] for k in (3,4,5))/max(len(P),1):.2f}  "
            f"origin R {ocount['R']/max(len(P),1):.2f}  "
            f"exploit depth1 {m(expl[1]):.2f} (n={len(expl[1])}) depth2+ {m(expl[2]):.2f} (n={len(expl[2])})  "
            f"copies-weighted depth1 {wd1:.2f}  "
            f"median age: parasite species {statistics.median(sages_p) if sages_p else float('nan'):.0f}, parasite lineages {statistics.median(ages_p) if ages_p else float('nan'):.0f}, replicator species {statistics.median(ages_r) if ages_r else float('nan'):.0f}")
    print(line, flush=True)
