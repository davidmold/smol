#!/usr/bin/env python3
"""Borrower partnerships: molecules that copy without a copy loop of their own.

    partners.py RUN_DIR [RUN_DIR ...]

A *borrower* is a species that makes exact copies of its partner (copies.tsv rows with
product == template) although its own sequence has no '=' (so it must be running its partner's code, e.g. via ^A). For each run:

  - first:   step the first borrower species is first seen copying
  - est:     first 10k window where borrower<->non-borrower copy events are >= 10% of all copying
  - origin:  the ancestry of the first established borrower, traced back through templates
             (species.tsv's template parent) to the nearest ancestor that has '='
  - share:   borrower<->non-borrower share of copy events per 100k window
"""
import sys, os
from collections import defaultdict

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


def ropen(p, mode="r", **k):
    """open p, or p.gz if p was compressed (archived runs)."""
    import gzip, os as _os
    if not _os.path.exists(p) and _os.path.exists(p + ".gz"):
        return gzip.open(p + ".gz", "rb" if "b" in mode else "rt", **k)
    return open(p, mode, **k)

def rexists(p):
    import os as _os
    return _os.path.exists(p) or _os.path.exists(p + ".gz")


def load_species(d):
    sp = {}
    for l in species_lines(os.path.join(d, "species.tsv")):
        p = l.rstrip("\n").split("\t")
        if p[0].isdigit():
            sp[int(p[0])] = (int(p[1]), int(p[2]), int(p[3]), p[-1])
    return sp

def analyse(d):
    sp = load_species(d)
    seq = lambda i: sp[i][3] if i in sp else ""
    rows = []
    for l in ropen(os.path.join(d, "copies.tsv")):
        p = l.split("\t")
        if p[0].isdigit() and p[2] == p[3]:  # exact copies of the template only
            rows.append((int(p[0]), int(p[1]), int(p[2]), int(p[4])))
    # borrowers = species that really made exact copies without a '=' of their own
    bset = {act for _, act, _, _ in rows if act in sp and "=" not in seq(act)}
    borrower = lambda i: i in bset
    win10 = defaultdict(lambda: [0, 0])      # step//10k -> [pair events, all events]
    win100 = defaultdict(lambda: [0, 0])
    first_copy = {}                           # borrower species -> first step seen copying
    last = 0
    for step, act, pas, n in rows:
        last = step
        ba, bp = borrower(act), borrower(pas)
        pair = ba != bp and (ba or bp)
        for w, k in ((win10, step // 10000), (win100, step // 100000)):
            w[k][1] += n
            if pair:
                w[k][0] += n
        if ba and act not in first_copy:
            first_copy[act] = step
    first = min(first_copy.values()) if first_copy else None
    est = next((k * 10000 for k in sorted(win10) if win10[k][1] > 0 and win10[k][0] / win10[k][1] >= 0.10), None)
    origin = []
    if est is not None:
        # the borrower that copied most in the establishing window's neighbourhood
        cand = [b for b, s in first_copy.items() if s <= est + 10000]
        b0 = min(cand, key=lambda b: first_copy[b]) if cand else None
        i, hops = b0, 0
        while i in sp and hops < 40:
            born, pa, pp, s = sp[i]
            origin.append((born, s))
            if "=" in s:
                break
            i, hops = pp, hops + 1
    shares = [(k * 100, round(100 * v[0] / v[1])) for k, v in sorted(win100.items()) if v[1]]
    alive = "EXTINCT" not in open(os.path.join(d, "log.txt")).read()
    return dict(run=d, alive=alive, last=last, first=first, est=est, origin=origin, shares=shares)

def mutation(a, b):
    """Describe how child b differs from parent a, roughly."""
    if len(a) == len(b):
        diff = [(x, y) for x, y in zip(a, b) if x != y]
        return "subst " + ",".join(f"{x}>{y}" for x, y in diff) if diff else "same"
    if len(b) < len(a):
        return f"shorter by {len(a) - len(b)}" + (" (truncated)" if a.startswith(b) or a.endswith(b) or b in a else "")
    return f"longer by {len(b) - len(a)}"

if __name__ == "__main__":
    for d in sys.argv[1:]:
        r = analyse(d.rstrip("/"))
        tag = "alive" if r["alive"] else "extinct"
        print(f"== {r['run']}  ({tag}, last {r['last']})")
        print(f"   first borrower copying: {r['first']}   established (>=10% of copying): {r['est']}")
        if r["origin"]:
            chain = list(reversed(r["origin"]))
            for k, (born, s) in enumerate(chain):
                how = "" if k == 0 else "   <- " + mutation(chain[k - 1][1], s)
                print(f"     {born:>8}  {s}{how}")
        print("   pair share by 100k:", " ".join(f"{s}%" for _, s in r["shares"]))
