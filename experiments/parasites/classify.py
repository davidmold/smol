#!/usr/bin/env python3
"""classify.py SRC_DIR STEP: classify the species present at STEP by what they did in the
following window (copies.tsv): R made exact copies (of itself or a partner); P never made a
copy but was copied exactly by someone (a parasite in the paper's sense); N neither.
Writes classes.tsv (species, class, molecules at STEP) and prints the molecule shares."""
import sys, gzip, collections
d, step = sys.argv[1], int(sys.argv[2])
mols = collections.Counter()
with gzip.open(f"{d}/checkpoint_{step}.txt.gz", "rt") as f:
    for l in f:
        if l.startswith("MOLS"): break
    for l in f:
        p = l.split()
        mols[int(p[3])] += 1
copier, copied = set(), set()
for l in open(f"{d}/copies.tsv"):
    p = l.split("\t")
    if len(p) < 5 or not p[0].isdigit(): continue
    a, pa, c = int(p[1]), int(p[2]), int(p[3])
    if c == pa or c == a:
        copier.add(a); copied.add(c)
cls = {}
for s in mols:
    cls[s] = "R" if s in copier else ("P" if s in copied else "N")
with open(f"{d}/classes.tsv", "w") as o:
    for s, n in mols.most_common():
        o.write(f"{s}\t{cls[s]}\t{n}\n")
tot = sum(mols.values()); by = collections.Counter()
for s, n in mols.items(): by[cls[s]] += n
print(d, tot, " ".join(f"{k}={by[k]} ({by[k]/tot:.0%})" for k in "RPN"),
      " species:", " ".join(f"{k}={sum(1 for s in cls if cls[s]==k)}" for k in "RPN"))
