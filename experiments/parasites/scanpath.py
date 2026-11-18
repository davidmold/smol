#!/usr/bin/env python3
"""scanpath.py RUNDIR... -> how did the first self-scanning replicator arise in each run?

For each run: find the first species (at >= MINN molecules in a popdy sample, up to 600k) that copies
itself exactly and is a self-scanner (copy ops per product symbol >= 1.5, as in compete.py). Then
walk its ancestry back through species.tsv (parent = whichever of the two recorded reactants is
nearer by edit distance) to the nearest ancestor that is a non-scanning self-copier. Every species on
the path in between is reported with: can it copy itself, scan or not, peak abundance, how long it
was present, and how many times it was copied (copies.tsv, any copier).

If the paths contain no non-self-copiers, self-scan arose in one step from a working replicator and
PARASITEFREE's blocking of non-viable intermediates cannot be what stopped it.
"""
import sys, os, collections, types
H = os.path.expanduser("~")
C = types.SimpleNamespace()
_src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "compete.py")).read()
_src = _src.split('if mode == "setup":')[0].replace("mode, OUT = sys.argv[1], sys.argv[2]", "")
exec(_src, C.__dict__)
MINN, MAXSTEP = 20, 600000

_cache = {}
def classify(s):
    if s in _cache: return _cache[s]
    res = C.parse(C.probe([s], ["--self"]))
    selfcop = 0 in res.get((0, 0), set())
    e = C.expo(s) if selfcop else 0.0
    _cache[s] = (selfcop, e); return _cache[s]

for d in sys.argv[1:]:
    seq, par = {}, {}
    for l in open(f"{d}/species.tsv"):
        p = l.rstrip("\n").split("\t")
        if p[0].isdigit() and len(p) >= 5: seq[int(p[0])] = p[-1]; par[int(p[0])] = (int(p[2]), int(p[3]))
    peak = collections.Counter(); first, last = {}, {}
    bystep = collections.defaultdict(dict)
    for l in open(f"{d}/popdy001.dat"):
        try: st, sp, n = (int(v) for v in l.split(","))
        except ValueError: continue
        if st > MAXSTEP: break
        bystep[st][sp] = n
        peak[sp] = max(peak[sp], n); first.setdefault(sp, st); last[sp] = st
    copied = collections.Counter()
    for l in open(f"{d}/copies.tsv"):
        p = l.split("\t")
        if len(p) < 5: continue
        if int(p[0]) > MAXSTEP: break
        if p[2] == p[3]: copied[int(p[2])] += int(p[4])
    target = None
    for st in sorted(bystep):
        for sp, n in sorted(bystep[st].items(), key=lambda kv: -kv[1]):
            if n < MINN: break
            if sp not in seq: continue
            sc, e = classify(seq[sp])
            if sc and e >= 1.5: target = (sp, st); break
        if target: break
    print(f"=== {d}")
    if not target: print("  no self-scanner up to", MAXSTEP); continue
    sp, st = target
    print(f"  first self-scanner: species {sp} at {st} (born {first.get(sp)})  {seq[sp]}")
    path, cur = [], sp
    while True:
        a, b = par[cur]
        cands = [x for x in (a, b) if x in seq]
        if not cands: break
        p = min(cands, key=lambda x: C.lev(seq[x], seq[cur]))
        sc, e = classify(seq[p])
        path.append((p, sc, e))
        if sc and e < 1.5: break
        if len(path) > 40: break
        cur = p
    nsc = 0
    for p, sc, e in path:
        kind = ("self-copier, scans" if e >= 1.5 else "self-copier, no scan") if sc else "CANNOT copy itself"
        nsc += (not sc)
        print(f"   <- {p:>6}  {kind:22} peak {peak[p]:>5}  present {first.get(p, '-')}-{last.get(p, '-')}  copied {copied[p]:>6}  {seq[p]}")
    print(f"  non-self-copying intermediates on the path: {nsc}")
