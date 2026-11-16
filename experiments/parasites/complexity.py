#!/usr/bin/env python3
"""complexity.py RUNS_DIR [--top 12] > table
For every run under RUNS_DIR/{none,P,cull}/sN and every kept checkpoint (each 100k), react
the top-K species pairwise in smol-react (mutation off) and report, for the most abundant
self-replicating species (the top replicator):
  len     its length
  q       number of '?' (label tests: the discrimination checks the paper describes)
  ops     number of executable opcodes in it ($ > ^ ? = % } +)
  rt      steps of its self-reaction
  sel     selectivity: share of the other top-K species it binds (P>0.05) that it does NOT copy
plus the run's population and the share of molecules in the top-K species that self-replicate.
Prints one line per run and step, then medians per mode and step."""
import sys, os, gzip, glob, re, subprocess, collections, statistics, tempfile

D = sys.argv[1]
K = int(sys.argv[sys.argv.index("--top") + 1]) if "--top" in sys.argv else 12
REACT = os.path.expanduser("~/smol-src/target/release/smol-react")
OPS = set("$>^?=%}+")

def snapshot(path):
    cnt = collections.Counter(); strands = collections.defaultdict(collections.Counter)
    with gzip.open(path, "rt") as f:
        for l in f:
            if l.startswith("MOLS"): break
        for l in f:
            p = l.split()
            sp = int(p[3]); cnt[sp] += 1
            if p[5] == "F":
                strands[sp][p[4].replace(".", "")] += 1
    seq = {sp: strands[sp].most_common(1)[0][0] for sp in strands}
    return cnt, seq

def react(seqs):
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        for i, s in enumerate(seqs):
            f.write(f"R m{i} {s}\n")
        name = f.name
    out = subprocess.run([REACT, name, "5000"], capture_output=True, text=True).stdout
    os.unlink(name)
    res = {}   # (x, y) -> (p, set of exact products, steps)
    cur = None
    for l in out.splitlines():
        m = re.match(r"(m\d+)\s+initiates with (m\d+)\s+P\(bind\)=([\d.]+).*steps=(\d+)", l)
        if m:
            cur = (int(m.group(1)[1:]), int(m.group(2)[1:])); res[cur] = [float(m.group(3)), set(), int(m.group(4))]
            continue
        m = re.match(r"\s+-> (?:product|passive made) (m\d+)", l)
        if m and cur:
            res[cur][1].add(int(m.group(1)[1:]))
    return res

rows = []
for d in sorted(glob.glob(f"{D}/*/s*")):
    mode, run = d.split("/")[-2], d.split("/")[-1]
    for ck in sorted(glob.glob(f"{d}/checkpoint_*00000.txt.gz"), key=lambda p: int(re.search(r"_(\d+)\.", p).group(1))):
        step = int(re.search(r"_(\d+)\.", ck).group(1))
        if step == 0: continue
        cnt, seq = snapshot(ck)
        top = [sp for sp, _ in cnt.most_common() if sp in seq][:K]
        if len(top) < 2: continue
        seqs = [seq[sp] for sp in top]
        res = react(seqs)
        selfrep = [i for i in range(len(top)) if i in res.get((i, i), [0, set()])[1]]
        pop = sum(cnt.values())
        srshare = sum(cnt[top[i]] for i in selfrep) / pop
        if not selfrep:
            rows.append((mode, run, step, pop, srshare, None)); continue
        i = selfrep[0]; s = seqs[i]
        bound = [j for j in range(len(top)) if j != i and max(res.get((i, j), [0])[0], res.get((j, i), [0])[0]) > 0.05]
        copied = [j for j in bound if j in res.get((i, j), [0, set()])[1] or j in res.get((j, i), [0, set()])[1]]
        sel = 1 - len(copied) / len(bound) if bound else float("nan")
        rows.append((mode, run, step, pop, srshare, (len(s), s.count("?"), sum(c in OPS for c in s), res[(i, i)][2], sel)))
        print(f"{mode:5} {run:5} {step:8d} pop {pop:6d} selfrep_share {srshare:.2f}  len {len(s):3d} q {s.count('?'):2d} ops {sum(c in OPS for c in s):3d} rt {res[(i,i)][2]:5d} sel {sel:.2f}  {s[:60]}", flush=True)

print("\n== medians per mode and step (runs with a self-replicating top-K species)")
print(f"{'mode':5} {'step':>8} {'n':>3} {'pop':>6} {'len':>5} {'q':>4} {'ops':>5} {'rt':>6} {'sel':>5}")
by = collections.defaultdict(list)
for r in rows:
    if r[5]: by[(r[0], r[2])].append(r)
for (mode, step) in sorted(by, key=lambda k: (k[0], k[1])):
    v = by[(mode, step)]
    med = lambda k: statistics.median([x[5][k] for x in v if x[5][k] == x[5][k]] or [float("nan")])
    print(f"{mode:5} {step:8d} {len(v):3d} {statistics.median([x[3] for x in v]):6.0f} {med(0):5.0f} {med(1):4.1f} {med(2):5.1f} {med(3):6.0f} {med(4):5.2f}")
