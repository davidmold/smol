#!/usr/bin/env python3
"""gm_parasites.py TS_DIR... -> what do gene-parasites look like? For each window (classes.tsv from
classify_gm.py), every gene-parasite species (P) is compared with the commonest working gene (G)
of the same window: length difference, and where the first and last differences fall (head =
before the coded part, i.e. the first 16 letters of the seed design; coded part after that), and
whether the gene's head (first 4 letters, the tag the machine's ?-test checks) and last 6 letters
are intact. Molecule-weighted shares over all windows."""
import sys, glob, os, collections, statistics as S
HEAD = 16
tot = collections.Counter(); n = 0; dl = []
for ts in sys.argv[1:]:
    for w in glob.glob(f"{ts}/s*/w*k"):
        if not os.path.exists(f"{w}/classes.tsv"): continue
        seq = {}
        for l in open(f"{w}/species.tsv"):
            p = l.rstrip("\n").split("\t")
            if p[0].isdigit() and len(p) >= 5: seq.setdefault(int(p[0]), p[-1])
        cl = [(int(a), b, int(c)) for a, b, c in (l.split() for l in open(f"{w}/classes.tsv"))]
        gs = [(c, a) for a, b, c in cl if b == "G" and a in seq]
        if not gs: continue
        g = seq[max(gs)[1]]
        for a, b, c in cl:
            if b != "P" or a not in seq: continue
            s = seq[a]; n += c
            dl.extend([len(s) - len(g)] * min(c, 50))
            i = next((k for k in range(min(len(s), len(g))) if s[k] != g[k]), min(len(s), len(g)))
            j = next((k for k in range(min(len(s), len(g))) if s[-1 - k] != g[-1 - k]), min(len(s), len(g)))
            tot["same length"] += c * (len(s) == len(g))
            tot["shorter"] += c * (len(s) < len(g)); tot["longer"] += c * (len(s) > len(g))
            tot["head tag intact (first 4)"] += c * (s[:4] == g[:4])
            tot["tail intact (last 6)"] += c * (s[-6:] == g[-6:])
            tot["first difference in head"] += c * (i < HEAD)
            tot["all differences in coded part"] += c * (i >= HEAD)
            tot["differs only in one stretch <= 4 letters (point-like)"] += c * (len(s) + len(g) - i - j <= 8 and abs(len(s) - len(g)) <= 2)
print(f"gene-parasite molecules: {n}; length difference from the commonest working gene: median {S.median(dl)}, "
      f"quartiles {S.quantiles(dl, n=4)}")
for k, v in tot.items(): print(f"  {k:55} {v / n:.2f}")
