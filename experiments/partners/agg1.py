#!/usr/bin/env python3
"""Aggregate partnership batch 1: one row per condition (pair x setup x mutation).

    agg1.py RUNS_DIR PAIRS_FILE

Per condition: runs alive at the end; median final shares of A, B, the copier, the parasite and
everything else; median fraction of the population that is A or B "exactly" (unmutated);
and over the last 20k steps, placed / attempted exact copies of A and B, and the share of A and
B molecules that died having had at least one placed exact copy made of them.
"""
import collections, glob, gzip, os, statistics, sys

def species_lines(p):
    """Lines of species.tsv and species.tsv.gz together (a run resumed from an edited checkpoint
    with GZIPLOGS has both: the copied list and the species born after the resume)."""
    import gzip as _g, os as _o
    for q in (p, p + ".gz"):
        if _o.path.exists(q):
            try:
                for l in (_g.open(q, "rt") if q.endswith(".gz") else open(q)):
                    yield l
            except EOFError:
                pass


runs, pairsf = sys.argv[1], sys.argv[2]
P = {}
for l in open(pairsf):
    if l.startswith("#") or not l.strip():
        continue
    p = l.split()
    P[p[0]] = {"A": p[1], "B": p[2], "par": p[3]}
COP = "RERERERE^BQQ$=?>%}"

def ropen(p):
    return open(p) if os.path.exists(p) else gzip.open(p + ".gz", "rt")

def one(d, names):
    seq = {}
    for l in species_lines(os.path.join(d, "species.tsv")):
        p = l.rstrip("\n").split("\t")
        if p[0].isdigit():
            seq[int(p[0])] = p[-1]
    inv = {v: k for k, v in names.items()}
    lab = lambda s: inv.get(seq.get(s, ""), "other")
    last = None; pop = collections.Counter()
    for l in ropen(os.path.join(d, "popdy001.dat")):
        t, s, c = l.split(","); t = int(t)
        if t != last:
            last = t; pop = collections.Counter()
        pop[lab(int(s))] += int(c)
    alive = "EXTINCT" not in open(os.path.join(d, "log.txt")).read()
    n = sum(pop.values())
    share = {k: pop[k] / n for k in pop} if n and alive else {}
    att = collections.Counter(); pl = collections.Counter()
    rows = [l.rstrip("\n").split("\t") for l in ropen(os.path.join(d, "copies.tsv"))]
    end = int(rows[-1][0]) if rows else 0
    for r in rows:
        if int(r[0]) > end - 20000 and r[2] == r[3]:
            k = lab(int(r[2])); att[k] += int(r[4]); pl[k] += int(r[6])
    died = collections.Counter(); cop = collections.Counter()
    lv = os.path.join(d, "lives.tsv")
    if os.path.exists(lv) or os.path.exists(lv + ".gz"):
        for l in ropen(lv):
            r = l.split("\t")
            if int(r[0]) > end - 20000:
                k = lab(int(r[1])); died[k] += int(r[2]); cop[k] += int(r[6])
    return alive, last, share, att, pl, died, cop

conds = collections.defaultdict(list)
for d in sorted(glob.glob(os.path.join(runs, "*", "seed*"))):
    conds[os.path.basename(os.path.dirname(d))].append(d)
print(f"{'condition':<18} {'alive':>5} {'A':>5} {'B':>5} {'cop':>5} {'par':>5} {'other':>5}  {'A placed/att':>12} {'B placed/att':>12} {'A dead copied':>13} {'B dead copied':>13}")
for c in sorted(conds):
    pair = c.split("_")[0]
    names = dict(P[pair]); names["cop"] = COP
    if "par" not in c:
        names.pop("par")
    res = [one(d, names) for d in conds[c]]
    alive = sum(r[0] for r in res)
    med = lambda k: statistics.median([r[2].get(k, 0) for r in res if r[0]]) if alive else 0
    def ratio(k, i, j):
        a = sum(r[i][k] for r in res); b = sum(r[j][k] for r in res)
        return f"{a/b:.3f}" if b else "-"
    print(f"{c:<18} {alive:>2}/{len(res):<2} {med('A'):>5.2f} {med('B'):>5.2f} {med('cop'):>5.2f} {med('par'):>5.2f} {med('other'):>5.2f}  {ratio('A',4,3):>12} {ratio('B',4,3):>12} {ratio('A',6,5):>13} {ratio('B',6,5):>13}")
