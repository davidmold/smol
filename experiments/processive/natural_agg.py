#!/usr/bin/env python3
"""Score the naturalistic introductions (natural.py; run on aida).

    natural_agg.py [--scan]     (--scan first runs transition_scan.py on finished runs)

Per run, from scan.json (transition_scan --every 20000 --min 2): established = processive share
>= 10% at the last sample (the endpoint test used for sizes.py, not the scan's 3-sample rule).
From trace.tsv (the newcomers and their descendants):
  founder's first encounter: none (decayed unbound) / P-resident (bound as the passive by a
  resident, i.e. a resident runs its code on the newcomer) / A-resident (the newcomer runs its
  code on a resident) / kin (bound to another tagged molecule)
  founder copied: an exact copy of the founder was placed; tagged product: any placed
  descendant (exact or mutant); lost: a product of a newcomer's reaction found no room;
  died bound: the founder died while in a reaction;
  lineage: alive at the end, or how many steps after introduction it went extinct;
  '>max' when the trace hit TRACEMAX (the lineage spread).
"""
import glob, json, os, statistics as st, subprocess, sys
from collections import Counter

H = os.path.expanduser("~/smol"); R = f"{H}/runs/natural"; CACHE = "/tmp/c_natural.json"

def step_done(d, tail):
    want = int(open(f"{d}/intro.txt").read().split()[2]) + 200000
    return f"Step {want} done" in tail

def trace(d):
    f = f"{d}/trace.tsv"
    if not os.path.exists(f):
        return None
    founders, alive, first, fate = [], set(), {}, {}
    nlines = 0; ext = None; maxalive = 0
    for l in open(f):
        nlines += 1
        p = l.rstrip("\n").split("\t"); t = int(p[0]); ev = p[1]
        if ev == "tag":
            k = (p[2], p[3]); founders.append(k); alive.add(k)
            fate[k] = dict(first=None, copied=0, tagged_prod=0, lost=0, died=None, died_bound=False, binds=0)
        elif ev == "bind":
            k = (p[2], p[3])
            if k in fate:
                fate[k]["binds"] += 1
                if fate[k]["first"] is None:
                    fate[k]["first"] = "kin" if p[9] == "1" else ("P-resident" if p[5] == "P" else "A-resident")
        elif ev == "prod":
            mk, pk = (p[2], p[3]), (p[6], p[7]); exact, tagged = p[11], p[12] == "1"
            for k in (mk, pk):
                if k in fate:
                    if (exact == "M" and k == mk) or (exact == "P" and k == pk):
                        fate[k]["copied"] += 1
                    if tagged:
                        fate[k]["tagged_prod"] += 1
            if tagged:
                alive.add((p[13], p[0])); maxalive = max(maxalive, len(alive))
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
    return dict(founders=[fate[k] for k in founders], ext=ext, alive=len(alive), nlines=nlines, maxalive=maxalive)

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

print(f"{'recipient':10} {'seq':4} {'N':>2} {'runs':>4} {'estab':>5} {'lin_alive':>9} {'ext_age_med':>11} "
      f"| founders: {'first: none/P-res/A-res/kin':27} {'copied':>6} {'desc':>4} {'lost':>4} {'died_bound':>10}")
for key in sorted(G):
    v = G[key]; rec, seq, n = key
    est = sum(1 for s, _ in v if s >= 0.10)
    trs = [t for _, t in v if t]
    alive = sum(1 for t in trs if t["alive"] > 0 or t["nlines"] >= 100000)
    ext_ages = [t["ext"] - t["t0"] for t in trs if t["ext"] is not None]
    F = [f for t in trs for f in t["founders"]]
    c = Counter(f["first"] or "none" for f in F)
    first = f"{c['none']}/{c['P-resident']}/{c['A-resident']}/{c['kin']}"
    cop = sum(1 for f in F if f["copied"] > 0); desc = sum(1 for f in F if f["tagged_prod"] > 0)
    lost = sum(1 for f in F if f["lost"] > 0); db = sum(1 for f in F if f["died_bound"])
    print(f"{rec:10} {seq:4} {n:>2} {len(v):>4} {est:>5} {alive:>9} {(str(int(st.median(ext_ages))) if ext_ages else '-'):>11} "
          f"| {len(F):>3} founders {first:>15} {cop:>6} {desc:>4} {lost:>4} {db:>10}")
