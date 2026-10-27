#!/usr/bin/env python3
"""Score the naturalistic introductions (natural.py; run on aida).

    natural_agg.py [--scan]     (--scan first runs transition_scan.py on finished runs)

Per run, from scan.json (transition_scan --every 20000 --min 2): established = processive share
>= 10% at the last sample (the endpoint test used for sizes.py, not the scan's 3-sample rule).
From trace.tsv (the newcomers and their descendants), per founder:
  first encounter: none / P-resident (a resident binds it as the passive and runs its own code
  on it) / A-resident (it runs its code on a resident) / kin (another tagged molecule);
  fate (exclusive): never bound (decayed free) / died in 1st reaction (died bound, no product,
  that reaction never ended) / reacted, no desc (reactions but no placed descendant) / had desc;
  altered: its own sequence was changed by a reaction (species id changed) before it died;
  copied: an exact copy of it (species id as at that bind) was placed; lost: a product of its
  reaction found no room.
Per run: lineage alive at the end (or the trace hit TRACEMAX, i.e. it spread), else how many
steps after introduction it went extinct.
"""
import glob, json, os, statistics as st, subprocess, sys
from collections import Counter

H = os.path.expanduser("~/smol"); R = f"{H}/runs/natural"; CACHE = "/tmp/c_natural.json"

def step_done(d, tail):
    want = int(open(f"{d}/intro.txt").read().split()[2]) + 200000
    return f"Step {want} done" in tail

def trace(d):
    """Per founder: first encounter, exclusive fate, whether residents altered its sequence."""
    f = f"{d}/trace.tsv"
    if not os.path.exists(f):
        return None
    founders, alive, fate = [], set(), {}
    nlines = 0; ext = None
    for l in open(f):
        nlines += 1
        p = l.rstrip("\n").split("\t"); t = int(p[0]); ev = p[1]
        if ev == "tag":
            k = (p[2], p[3]); founders.append(k); alive.add(k)
            fate[k] = dict(sp0=p[4], first=None, altered=False, copied=0, desc=0, lost=0, ends=0,
                           prods=0, died=None, died_bound=False)
            continue
        if ev in ("bind", "end", "die"):
            k = (p[2], p[3])
            if k in fate and p[4] != fate[k]["sp0"]:
                fate[k]["altered"] = True
        if ev == "bind":
            k = (p[2], p[3])
            if k in fate and fate[k]["first"] is None:
                fate[k]["first"] = "kin" if p[9] == "1" else ("P-resident" if p[5] == "P" else "A-resident")
        elif ev == "end":
            k = (p[2], p[3])
            if k in fate:
                fate[k]["ends"] += 1
        elif ev == "prod":
            mk, pk = (p[2], p[3]), (p[6], p[7]); msp, psp, sp, tagged = p[4], p[8], p[10], p[12] == "1"
            for k, ksp in ((mk, msp), (pk, psp)):
                if k in fate:
                    fate[k]["prods"] += 1
                    if sp == ksp:  # an exact copy of the founder as it was in this reaction
                        fate[k]["copied"] += 1
                    if tagged:
                        fate[k]["desc"] += 1
            if tagged:
                alive.add((p[13], p[0]))
        elif ev == "lost":
            for k in ((p[2], p[3]), (p[6], p[7])):
                if k in fate:
                    fate[k]["lost"] += 1
        elif ev == "die":
            k = (p[2], p[3])
            if k in fate:
                fate[k]["died"] = t; fate[k]["died_bound"] = p[5] == "1"
            alive.discard(k)
            if not alive and ext is None:
                ext = t
    for k, v in fate.items():
        if v["first"] is None:
            v["fate"] = "never bound"
        elif v["died"] is not None and v["died_bound"] and v["ends"] == 0 and v["prods"] == 0:
            v["fate"] = "died in 1st reaction"
        elif v["desc"] == 0:
            v["fate"] = "reacted, no desc"
        else:
            v["fate"] = "had desc"
    return dict(founders=[fate[k] for k in founders], ext=ext, alive=len(alive), nlines=nlines)

runs = sorted(os.path.dirname(p) for p in glob.glob(f"{R}/*/*/n*/rep*/log.txt"))
if "--scan" in sys.argv:
    if not os.path.exists(CACHE) and os.path.exists("/tmp/cs.json"):
        subprocess.run(["cp", "/tmp/cs.json", CACHE])
    for d in runs:
        st0 = int(open(f"{d}/intro.txt").read().split()[2])
        t = open(f"{d}/log.txt", errors="replace").read()[-4000:]
        if (step_done(d, t) or "EXTINCT" in t) and not os.path.exists(f"{d}/scan.json"):
            subprocess.run([sys.executable, "-I", f"{H}/tools/transition_scan.py", d, "--every", "20000",
                            "--start", str(st0), "--min", "2", "--cache", CACHE, "--json", f"{d}/scan.json"],
                           stdout=subprocess.DEVNULL, check=True)

G = {}
for d in runs:
    if not os.path.exists(f"{d}/scan.json"):
        continue
    rec, seq, n, rep = d.split("/")[-4:]
    r = json.load(open(f"{d}/scan.json"))["rows"]
    share = r[-1]["proc_share"] if r else 0.0
    tr = trace(d)
    if tr:
        tr["t0"] = int(open(f"{d}/intro.txt").read().split()[2])
    G.setdefault((rec, seq, int(n[1:])), []).append((share, tr))

FATES = ["never bound", "died in 1st reaction", "reacted, no desc", "had desc"]
print(f"{'recipient':10} {'seq':4} {'N':>2} {'runs':>4} {'estab':>5} {'lin_alive':>9} {'ext_age_med':>11} | "
      f"{'founders':>8} {'first none/P/A/kin':>18} {'fate never/trap/nodesc/desc':>27} {'altered':>7} {'copied':>6} {'lost':>4}")
for key in sorted(G):
    v = G[key]; rec, seq, n = key
    est = sum(1 for s_, _ in v if s_ >= 0.10)
    trs = [t for _, t in v if t]
    alive = sum(1 for t in trs if t["alive"] > 0 or t["nlines"] >= 100000)
    ext_ages = [t["ext"] - t["t0"] for t in trs if t["ext"] is not None]
    F = [f for t in trs for f in t["founders"]]
    c = Counter(f["first"] or "none" for f in F); fc = Counter(f["fate"] for f in F)
    first = f"{c['none']}/{c['P-resident']}/{c['A-resident']}/{c['kin']}"
    fates = "/".join(str(fc[x]) for x in FATES)
    print(f"{rec:10} {seq:4} {n:>2} {len(v):>4} {est:>5} {alive:>9} {(str(int(st.median(ext_ages))) if ext_ages else '-'):>11} | "
          f"{len(F):>8} {first:>18} {fates:>27} {sum(f['altered'] for f in F):>7} "
          f"{sum(1 for f in F if f['copied'] > 0):>6} {sum(1 for f in F if f['lost'] > 0):>4}")
