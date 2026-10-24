#!/usr/bin/env python3
"""Aggregate lineage.py JSON outputs per condition: lineage_agg.py LIN_DIR

Per run, the organisation at 200k is classed from placed exact copies in the last 20k steps:
  partnership  A<->B cross-lineage copies >= 50% of all, each direction >= 5%
  autonomous   same-species self-copies >= 50%
  mixed        anything else (e.g. copier lineage or within-lineage copying dominant)
Also: share of each lineage's molecules (at >= 1% species) that are in species which made no
placed copies themselves ("passengers": copied, never copying), and the effective number of
species pairs carrying the B->A and A->B copies.
"""
import collections, glob, json, os, statistics, sys

d = sys.argv[1]
by = collections.defaultdict(list)
for f in sorted(glob.glob(f"{d}/*.json")):
    o = json.load(open(f))
    by[os.path.basename(f).rsplit("_seed", 1)[0]].append(o)
med = lambda xs: statistics.median(xs) if xs else float("nan")
print(f"{'condition':<12} {'alive':>5} {'A':>5} {'B':>5} {'cop':>5} {'par':>5} | {'A<->B':>5} {'self':>5} {'withinL':>7} | {'class P/Au/M':>12} | {'B pass':>6} {'A pass':>6} | {'effB->A':>7} {'effA->B':>7}")
for c, os_ in sorted(by.items()):
    al = [o for o in os_ if o["alive"]]
    def sh(o, k): return o["pop"].get(k, 0)
    def cp(o, pred): return sum(v for k, v in o["copies"].items() if pred(k))
    ab = [cp(o, lambda k: k in ("A->B", "B->A")) for o in al]
    se = [cp(o, lambda k: k == "self") for o in al]
    wl = [cp(o, lambda k: k in ("A->A", "B->B", "cop->cop", "par->par")) for o in al]
    cls = collections.Counter()
    for o in al:
        x = cp(o, lambda k: k in ("A->B", "B->A")); s = o["copies"].get("self", 0)
        if x >= .5 and o["copies"].get("A->B", 0) >= .05 and o["copies"].get("B->A", 0) >= .05:
            cls["P"] += 1
        elif s >= .5:
            cls["Au"] += 1
        else:
            cls["M"] += 1
    def passengers(o, lin):
        tot = sum(r["share"] for r in o["roles"] if r["lin"] == lin)
        pas = sum(r["share"] for r in o["roles"] if r["lin"] == lin and r["made"] == 0)
        return pas / tot if tot else float("nan")
    bp = [passengers(o, "B") for o in al]; apas = [passengers(o, "A") for o in al]
    eb = [o["cross"]["B->A"]["eff_pairs"] for o in al if o["cross"]["B->A"]["copies"]]
    ea = [o["cross"]["A->B"]["eff_pairs"] for o in al if o["cross"]["A->B"]["copies"]]
    nn = lambda xs: [x for x in xs if x == x]
    print(f"{c:<12} {len(al):>2}/{len(os_):<2} {med([sh(o,'A') for o in al]):>5.2f} {med([sh(o,'B') for o in al]):>5.2f} {med([sh(o,'cop') for o in al]):>5.2f} {med([sh(o,'par') for o in al]):>5.2f} | {med(ab):>5.2f} {med(se):>5.2f} {med(wl):>7.2f} | {cls['P']:>3}/{cls['Au']:>2}/{cls['M']:>2}     | {med(nn(bp)):>6.2f} {med(nn(apas)):>6.2f} | {med(eb):>7.1f} {med(ea):>7.1f}")
