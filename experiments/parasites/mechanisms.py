#!/usr/bin/env python3
"""mechanisms.py PARASITES_DIR > table

Fingerprints of tolerance and of "evolution evolving", at 1M and 2M in each published run
(the converted worlds and classes from prep.sh / classify.py, and the swap control runs).

For the top replicator S (most abundant self-copier) and the world's 10 commonest parasites
(class P), in the reaction probe (mutation off):
  tS        steps per copy when S copies itself
  ratio     median over parasites of (steps per parasite copy) / tS: below 1 means parasites are
            copied faster than kin; near 1 means their speed advantage is gone
  cps       copies of itself S makes per self-binding
  twoway    binding-weighted share of S-parasite encounters in which S is also copied
  expo      copy-operator executions per product symbol in S's self-reaction (1 = plain copy;
            about 2 = self-scan, which roughly doubles exposure to point mutation)
From the world snapshot:
  clusR     share of a replicator's occupied neighbours that are replicators, over the global
            replicator share (1 = random mixing; above 1 = replicators sit in patches)
From the 100k-step control runs that follow each snapshot (swap self1 / self2), copies.tsv:
  big       share of new products (not an exact copy) that differ from both reactants by more than
            2 edits: rearrangements, not point mutations
"""
import sys, os, re, gzip, glob, collections, subprocess, tempfile, statistics

P = sys.argv[1]
REACT = os.path.expanduser("~/smol-src/target/release/smol-react")

def lev(a, b, cap=3):
    if abs(len(a) - len(b)) > cap: return cap + 1
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i] + [0] * len(b)
        for j, cb in enumerate(b, 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb))
        if min(cur) > cap: return cap + 1
        prev = cur
    return prev[-1]

def species(d):
    s = {}
    for l in open(f"{d}/species.tsv"):
        p = l.rstrip("\n").split("\t")
        if p[0].isdigit(): s.setdefault(int(p[0]), p[-1])
    return s

def probe(seqs, flags):
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        for i, s in enumerate(seqs): f.write(f"R m{i} {s}\n")
        name = f.name
    out = subprocess.run([REACT, name, "5000"] + flags, capture_output=True, text=True).stdout
    os.unlink(name)
    return out

def parse(out):
    res = {}; cur = None
    for l in out.splitlines():
        m = re.match(r"(m\d+)\s+initiates with (m\d+)\s+P\(bind\)=([\d.]+).*steps=(\d+)", l)
        if m:
            cur = (int(m.group(1)[1:]), int(m.group(2)[1:]))
            res[cur] = dict(p=float(m.group(3)), steps=int(m.group(4)), prods=collections.Counter())
            continue
        m = re.match(r"\s+-> (?:product|passive made) (m\d+)", l)
        if m and cur: res[cur]["prods"][int(m.group(1)[1:])] += 1
    return res

def world(src, step):
    cls = {int(l.split()[0]): l.split()[1] for l in open(f"{src}/classes.tsv")}
    cnt = collections.Counter(); cell = {}; gx = gy = None
    with gzip.open(f"{src}/checkpoint_{step}.txt.gz", "rt") as f:
        for l in f:
            if l.startswith("GRIDX"): gx = int(l.split()[1])
            if l.startswith("GRIDY"): gy = int(l.split()[1])
            if l.startswith("MOLS"): break
        for l in f:
            p = l.split(); sp = int(p[3]); cnt[sp] += 1; cell[(int(p[1]), int(p[2]))] = sp
    return cls, cnt, cell, gx, gy

def clustering(cls, cell, gx, gy):
    isR = {c: cls.get(sp) == "R" for c, sp in cell.items()}
    glob_r = sum(isR.values()) / len(isR)
    fr = []
    for (x, y), r in isR.items():
        if not r: continue
        nb = [isR[((x + dx) % gx, (y + dy) % gy)] for dx in (-1, 0, 1) for dy in (-1, 0, 1)
              if (dx or dy) and ((x + dx) % gx, (y + dy) % gy) in isR]
        if nb: fr.append(sum(nb) / len(nb))
    return statistics.mean(fr) / glob_r if fr and glob_r else float("nan")

def rearrangements(rundirs):
    big = tot = 0
    for d in rundirs:
        seq = species(d)
        cache = {}
        for l in open(f"{d}/copies.tsv"):
            p = l.split("\t")
            if len(p) < 5 or not p[0].isdigit(): continue
            a, pa, c, n = int(p[1]), int(p[2]), int(p[3]), int(p[4])
            if c in (a, pa) or c not in seq: continue
            k = (a, pa, c)
            if k not in cache:
                sc = seq[c]
                cache[k] = min(lev(sc, seq.get(a, "")), lev(sc, seq.get(pa, "")))
            tot += n; big += n * (cache[k] > 2)
    return big / tot if tot else float("nan")

rows = collections.defaultdict(dict)
for step, ctl in ((1000000, "self1"), (2000000, "self2")):
    for src in sorted(glob.glob(f"{P}/s{step}/r*/src")):
        r = src.split("/")[-2]
        cls, cnt, cell, gx, gy = world(src, step)
        seq = species(src)
        top = [sp for sp, _ in cnt.most_common(12) if sp in seq]
        sr = parse(probe([seq[s] for s in top], ["--self"]))
        S = next((top[i] for i in range(len(top)) if sr.get((i, i), {}).get("prods", {}).get(i)), None)
        out = dict(clus=clustering(cls, cell, gx, gy),
                   big=rearrangements(glob.glob(f"{P}/swap/{r}/{ctl}/s*")))
        if S is not None:
            par = [sp for sp, _ in cnt.most_common() if cls.get(sp) == "P" and sp in seq][:10]
            res = parse(probe([seq[S]] + [seq[p] for p in par], ["--hub"]))
            selfr = res.get((0, 0)); cps = selfr["prods"][0] if selfr else 0
            tS = selfr["steps"] / cps if cps else float("nan")
            ratios = []; tw_num = tw_den = 0.0
            for j in range(1, len(par) + 1):
                e = res.get((0, j))
                if not e or e["p"] <= 0: continue
                if e["prods"][j]: ratios.append(e["steps"] / e["prods"][j] / tS)
                if e["prods"][j]:
                    tw_den += e["p"]; tw_num += e["p"] * (e["prods"][0] > 0)
            tr = probe([seq[S]], ["--steptrace", "--self"])
            ncopy = sum(1 for l in tr.splitlines() if l.startswith("T ") and l.split()[2] == "=")
            plen = sum(len(l.split()[2]) for l in tr.splitlines() if l.startswith("PROD"))
            out.update(len=len(seq[S]), tS=tS, ratio=statistics.median(ratios) if ratios else float("nan"),
                       cps=cps, twoway=tw_num / tw_den if tw_den else float("nan"),
                       expo=ncopy / plen if plen else float("nan"))
        rows[r][step] = out
        f = lambda k: out.get(k, float("nan"))
        print(f"{r:4} {step:8d} len {f('len'):4.0f} tS {f('tS'):5.1f} ratio {f('ratio'):4.2f} cps {f('cps'):3.0f} twoway {f('twoway'):4.2f} expo {f('expo'):4.2f} clusR {f('clus'):4.2f} big {f('big'):.3f}", flush=True)

print("\n== medians, 1M -> 2M (runs with values at both)")
for k in ("tS", "ratio", "cps", "twoway", "expo", "clus", "big"):
    a = [rows[r][1000000].get(k) for r in rows if 2000000 in rows[r] and 1000000 in rows[r]]
    b = [rows[r][2000000].get(k) for r in rows if 2000000 in rows[r] and 1000000 in rows[r]]
    pr = [(x, y) for x, y in zip(a, b) if x is not None and y is not None and x == x and y == y]
    if pr:
        print(f"  {k:7} {statistics.median([x for x,_ in pr]):7.3f} -> {statistics.median([y for _,y in pr]):7.3f}   up in {sum(y>x for x,y in pr)}/{len(pr)} runs")
