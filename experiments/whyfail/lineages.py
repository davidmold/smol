#!/usr/bin/env python3
"""Processive lineages in one run, from its saved histories (run on aida).

    lineages.py RUN_DIR T0 T_END OUT.json [--cache C]

1. popdy: every species that ever had >= 3 molecules (candidates), with its count every 100 steps.
2. species.tsv(.gz): sequence, birth step, parents (pa = maker, pp = template) of candidates,
   and of their parents.
3. probe candidates (smol-react --self, mutation off, 5,000 steps; cached): processive = binds
   itself and makes >= 3 exact self-copies in one reaction (the transition_scan definition).
4. lineages: processive species linked to a processive parent belong to that parent's lineage;
   a processive species whose parents are not processive founds a new lineage (an independent
   origin). Species that never reached 3 molecules are not followed.
Per lineage: founder (id, sequence, birth, its parents' sequences and whether they were
processive), first step with >= 1/5/20 molecules, peak, share at each 10k, member species
(with birth and peak), fate (established = >= 10% of molecules for >= 3 consecutive 10k
samples; else the last step it had molecules).
"""
import gzip, json, os, re, subprocess, sys, tempfile
from collections import defaultdict

run, T0, TEND, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
cache_f = sys.argv[sys.argv.index("--cache") + 1] if "--cache" in sys.argv else None
REACT = os.path.expanduser("~/smol/target/release/smol-react")

def lines(p):
    for q in (p, p + ".gz"):
        if os.path.exists(q):
            try:
                for l in (gzip.open(q, "rt") if q.endswith(".gz") else open(q)):
                    yield l
            except EOFError:
                pass

# 1. popdy
counts = defaultdict(dict); total = defaultdict(int); maxc = defaultdict(int)
for l in lines(f"{run}/popdy001.dat"):
    t, s, c = l.split(","); t = int(t)
    if t < T0 or t > TEND:
        continue
    s = int(s); c = int(c); total[t] += c
    if c > maxc[s]:
        maxc[s] = c
    if c >= 1:
        counts[s][t] = c  # kept for all; trimmed below
cand = {s for s, m in maxc.items() if m >= 3}
counts = {s: counts[s] for s in cand}

# 2. species table for candidates (+ parents)
info = {}
for l in lines(f"{run}/species.tsv"):
    i = l.find("\t")
    if i > 0 and l[:i].isdigit() and int(l[:i]) in cand:
        p = l.rstrip("\n").split("\t")
        info[int(p[0])] = dict(born=int(p[1]), pa=int(p[2]), pp=int(p[3]), seq=p[4])
par = {v["pa"] for v in info.values()} | {v["pp"] for v in info.values()}
pseq = {}
for l in lines(f"{run}/species.tsv"):
    i = l.find("\t")
    if i > 0 and l[:i].isdigit() and int(l[:i]) in par:
        p = l.rstrip("\n").split("\t"); pseq[int(p[0])] = p[4]

# 3. probe
cache = json.load(open(cache_f)) if cache_f and os.path.exists(cache_f) else {}
todo = sorted({v["seq"] for v in info.values()} | set(pseq.values()) - set(cache))
pat = re.compile(r"(\S+)\s+initiates with (\S+)\s+P\(bind\)=([\d.]+).*steps=(\d+)")
for i in range(0, len(todo), 400):
    chunk = todo[i:i + 400]
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        for k, sq in enumerate(chunk):
            f.write(f"R m{k} {sq}\n")
    r = subprocess.run([REACT, f.name, "5000", "--self"], capture_output=True, text=True)
    os.unlink(f.name)
    if r.returncode != 0:
        sys.exit(f"smol-react failed: {r.stderr[-300:]}")
    res = {}; cur = None
    for l in r.stdout.split("\n"):
        m = pat.match(l)
        if m:
            cur = m.group(1); res[cur] = [float(m.group(3)), int(m.group(4)), 0]
        elif "-> product" in l and cur and l.split()[2] == cur:
            res[cur][2] += 1
    for k, sq in enumerate(chunk):
        if f"m{k}" not in res:
            sys.exit("missing probe result")
        cache[sq] = res[f"m{k}"]
if cache_f:  # atomic: a reader never sees a half-written cache
    json.dump(cache, open(cache_f + ".tmp", "w")); os.replace(cache_f + ".tmp", cache_f)
proc = lambda sq: sq in cache and cache[sq][0] > 0 and cache[sq][2] >= 3

# 4. lineages
P = {s for s in info if proc(info[s]["seq"])}
root = {}
def find(s, depth=0):
    if s in root:
        return root[s]
    v = info[s]; r = s
    for q in (v["pp"], v["pa"]):  # template first: a copy descends from what was copied
        if q in P and q != s and depth < 500 and info[q]["born"] <= v["born"]:
            r = find(q, depth + 1); break
    root[s] = r; return r
for s in sorted(P, key=lambda s: info[s]["born"]):
    find(s)
lin = defaultdict(list)
for s in P:
    lin[root[s]].append(s)
steps = sorted(total)
samples = [t for t in steps if (t - T0) % 10000 == 0 and t > T0]
L = []
for r, mem in lin.items():
    n = defaultdict(int)
    for s in mem:
        for t, c in counts[s].items():
            n[t] += c
    first = lambda k: next((t for t in steps if n.get(t, 0) >= k), None)
    share = {t: n.get(t, 0) / total[t] for t in samples}
    run_len = 0; est = None
    for i, t in enumerate(samples):
        run_len = run_len + 1 if share[t] >= 0.10 else 0
        if est is None and run_len >= 3:
            est = samples[i - 2]
    last = max((t for t in steps if n.get(t, 0) > 0), default=None)
    v = info[r]
    L.append(dict(root=r, seq=v["seq"], len=len(v["seq"]), born=v["born"], probe=cache[v["seq"]],
                  maker=v["pa"], maker_seq=pseq.get(v["pa"]), maker_proc=proc(pseq.get(v["pa"], "")),
                  template=v["pp"], template_seq=pseq.get(v["pp"]), template_proc=proc(pseq.get(v["pp"], "")),
                  first1=first(1), first5=first(5), first20=first(20), peak=max(n.values(), default=0),
                  peak_step=max(n, key=n.get) if n else None, est=est, last=last,
                  share={str(t): round(share[t], 4) for t in samples},
                  members=sorted([dict(id=s, seq=info[s]["seq"], born=info[s]["born"],
                                       peak=max(counts[s].values())) for s in mem], key=lambda m: -m["peak"])[:20],
                  n_members=len(mem), member_ids=sorted(mem)))
L.sort(key=lambda x: -x["peak"])
json.dump(dict(run=run, T0=T0, TEND=TEND, candidates=len(cand), processive=len(P), lineages=L), open(out, "w"), indent=1)
print(f"{run}: {len(cand)} candidates, {len(P)} processive, {len(L)} lineages; "
      + "; ".join(f"root {x['root']} len {x['len']} born {x['born']} peak {x['peak']} est {x['est']} last {x['last']}" for x in L[:5]))
