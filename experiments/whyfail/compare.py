#!/usr/bin/env python3
"""What first differs between processive lineages that establish and those that disappear?
(run on aida after lineages.py)

    compare.py LINEAGE_JSON... > compare.txt   (writes compare_rows.json too)

For every lineage that reached >= 5 molecules (a foothold), in the window W = [t5, t5 + 5000]
where t5 is its first step with >= 5 molecules:
  genotype     founder length; probe: P(bind) with itself, exact self-copies per reaction
  copying      copies.tsv rows in W where a lineage member is involved, by kind:
                 kin   member (active) on member (passive)
                 used  resident runs its code on a member (member passive)
                 uses  member runs its code on a resident
               products of each kind that are lineage members (prod in the lineage), and placed
  space        grid.snap frames in W: of the 8 neighbours of member cells, fraction empty / member /
               resident; fraction of members bound; members' largest connected cluster
  growth       members in the frames at t5, t5+2000, t5+5000
  later mutation (established lineages): the member most abundant at establishment, its birth
               and length relative to the founder
Established = lineages.py's rule (>= 10% of molecules for >= 3 consecutive 10k samples).
"""
import gzip, json, os, statistics as st, struct, sys, zlib
from collections import Counter

W = 5000

def frames(run, t0, t1):
    f = open(f"{run}/grid.snap", "rb")
    while True:
        h = f.read(16)
        if len(h) < 16:
            return
        s, gx, gy, L = struct.unpack("<4I", h)
        if s < t0:
            f.seek(L, 1); continue
        if s > t1:
            return
        raw = zlib.decompress(f.read(L))
        yield s, gx, gy, struct.unpack(f"<{gx * gy}I", raw)

def copies(run, t0, t1):
    for p in (f"{run}/copies.tsv", f"{run}/copies.tsv.gz"):
        if os.path.exists(p):
            try:
                for l in (gzip.open(p, "rt") if p.endswith(".gz") else open(p)):
                    r = l.split("\t"); t = int(r[0])
                    if t < t0:
                        continue
                    if t > t1:
                        return
                    yield int(r[1]), int(r[2]), int(r[3]), int(r[4]), int(r[6]) if len(r) > 6 else 0
            except EOFError:
                return

rows = []
for jf in sys.argv[1:]:
    J = json.load(open(jf)); run = J["run"]
    for x in J["lineages"]:
        if x["first5"] is None:
            continue
        M = set(x["member_ids"]); t5 = x["first5"]
        k = Counter()
        for a, p, pr, att, placed in copies(run, t5, t5 + W):
            ia, ip = a in M, p in M
            if not (ia or ip):
                continue
            kind = "kin" if ia and ip else ("used" if ip else "uses")
            k[kind] += att; k[kind + "_placed"] += placed
            if pr in M:
                k[kind + "_member"] += att; k[kind + "_member_placed"] += placed
        nb = Counter(); bound = n_mem = 0; sizes = {}; clus = []
        for s, gx, gy, g in frames(run, t5 - 499, t5 + W):
            cells = [i for i, v in enumerate(g) if v and (v & 0x7fffffff) in M]
            sizes[s] = len(cells)
            n_mem += len(cells); bound += sum(1 for i in cells if g[i] & 0x80000000)
            cs = set(cells); seen = set(); best = 0
            for i in cells:
                y, xx = divmod(i, gx)
                for dy in (-1, 0, 1):
                    for dx in (-1, 0, 1):
                        if dx or dy:
                            j = ((y + dy) % gy) * gx + (xx + dx) % gx
                            v = g[j]
                            nb["empty" if not v else ("member" if (v & 0x7fffffff) in M else "resident")] += 1
                if i not in seen:  # cluster size by flood fill (8-neighbourhood)
                    stack = [i]; seen.add(i); n = 0
                    while stack:
                        c = stack.pop(); n += 1; cy, cx = divmod(c, gx)
                        for dy in (-1, 0, 1):
                            for dx in (-1, 0, 1):
                                j = ((cy + dy) % gy) * gx + (cx + dx) % gx
                                if j in cs and j not in seen:
                                    seen.add(j); stack.append(j)
                    best = max(best, n)
            if cells:
                clus.append(best / len(cells))
        tn = sum(nb.values()) or 1
        fs = sorted(sizes)
        at = lambda dt: next((sizes[s] for s in fs if s >= t5 + dt), None)
        dom = None
        if x["est"]:
            # member most abundant at establishment: approximate by the member with the highest peak
            m = x["members"][0]; dom = dict(len=len(m["seq"]), born_after_root=m["born"] - x["born"],
                                            born_after_t5=m["born"] - t5, is_root=m["id"] == x["root"])
        rows.append(dict(run=run.split("/runs/")[-1], root=x["root"], est=x["est"] is not None, len=x["len"],
                         pbind=x["probe"][0], selfcopies=x["probe"][2], t5=t5, peak=x["peak"],
                         copying=dict(k), empty=nb["empty"] / tn, member_nb=nb["member"] / tn,
                         resident_nb=nb["resident"] / tn, bound=bound / max(n_mem, 1),
                         cluster=st.mean(clus) if clus else None,
                         n0=at(0), n2k=at(2000), n5k=at(5000), dominant=dom))

json.dump(rows, open("compare_rows.json", "w"), indent=1)
S = [r for r in rows if r["est"]]; F = [r for r in rows if not r["est"]]
print(f"lineages with >= 5 molecules: {len(rows)} ({len(S)} established, {len(F)} not)\n")
def med(v):
    v = [x for x in v if x is not None]
    return f"{st.median(v):.3g}" if v else "-"
def frac(r, a, b):
    c = r["copying"]; return c.get(a, 0) / c[b] if c.get(b) else None
feat = [("founder length", lambda r: r["len"]), ("probe P(bind)", lambda r: r["pbind"]),
        ("probe self-copies", lambda r: r["selfcopies"]),
        ("members at t5", lambda r: r["n0"]), ("members at t5+2k", lambda r: r["n2k"]), ("members at t5+5k", lambda r: r["n5k"]),
        ("copy attempts: kin", lambda r: r["copying"].get("kin", 0)),
        ("copy attempts: used by residents", lambda r: r["copying"].get("used", 0)),
        ("copy attempts: using residents", lambda r: r["copying"].get("uses", 0)),
        ("kin share of member reactions", lambda r: (r["copying"].get("kin", 0) / max(1, sum(r["copying"].get(q, 0) for q in ("kin", "used", "uses"))))),
        ("used: product is a member", lambda r: frac(r, "used_member", "used")),
        ("kin: product placed", lambda r: frac(r, "kin_placed", "kin")),
        ("neighbours empty", lambda r: r["empty"]), ("neighbours member", lambda r: r["member_nb"]),
        ("neighbours resident", lambda r: r["resident_nb"]), ("members bound", lambda r: r["bound"]),
        ("largest cluster share", lambda r: r["cluster"])]
print(f"{'feature (window t5..t5+5k)':38} {'established':>12} {'failed':>12}")
for name, f in feat:
    print(f"{name:38} {med(f(r) for r in S):>12} {med(f(r) for r in F):>12}")
print("\nestablished lineages:")
for r in S:
    print(f"  {r['run']:28} len {r['len']:3} self-copies {r['selfcopies']:3} n0/2k/5k {r['n0']}/{r['n2k']}/{r['n5k']} "
          f"empty {r['empty']:.2f} member {r['member_nb']:.2f} dominant {r['dominant']}")
