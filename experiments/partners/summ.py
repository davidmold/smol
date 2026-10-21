#!/usr/bin/env python3
"""Summarise a partnership run: summ.py RUN_DIR [NAME=SEQ ...]

Shares of named sequences at the last popdy step, placed vs attempted copies made of each
(template = product) over the last 20k steps, and from lives.tsv the share of molecules that
died having been copied at least once (successful reproducers).
"""
import sys, collections, os
d = sys.argv[1]
named = dict(kv.split("=", 1) for kv in sys.argv[2:])
seq = {}
for l in open(os.path.join(d, "species.tsv")):
    p = l.rstrip("\n").split("\t")
    if p[0].isdigit():
        seq[int(p[0])] = p[-1]
name_of = {v: k for k, v in named.items()}
lab = lambda s: name_of.get(seq.get(s, ""), "other")
last = None; pop = collections.Counter()
for l in open(os.path.join(d, "popdy001.dat")):
    t, s, c = map(int, l.split(","))
    if t != last:
        last = t; pop = collections.Counter()
    pop[lab(s)] += c
n = sum(pop.values())
out = [f"{os.path.basename(d.rstrip('/'))} step {last} n={n}"]
out.append(" share " + " ".join(f"{k}:{v/n:.2f}" for k, v in pop.most_common()) if n else " extinct")
rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(d, "copies.tsv"))]
if rows:
    end = int(rows[-1][0])
    att = collections.Counter(); pl = collections.Counter()
    for r in rows:
        if int(r[0]) > end - 20000 and r[2] == r[3]:
            att[lab(int(r[2]))] += int(r[4]); pl[lab(int(r[2]))] += int(r[6]) if len(r) > 6 else 0
    out.append(" exact copies of (last 20k) " + " ".join(f"{k}:{pl[k]}/{att[k]}" for k in att))
L = collections.defaultdict(lambda: [0, 0])
if os.path.exists(os.path.join(d, "lives.tsv")):
    for l in open(os.path.join(d, "lives.tsv")):
        r = l.split("\t")
        L[lab(int(r[1]))][0] += int(r[2]); L[lab(int(r[1]))][1] += int(r[6])
    out.append(" died-having-been-copied " + " ".join(f"{k}:{v[1]}/{v[0]}" for k, v in L.items() if v[0]))
print("\n".join(out))
