#!/usr/bin/env python3
"""gm_scan2.py -> where in the machine can substitutions land, and what do indels do?
(Uses gm_scan.py's react() and seed.) For each machine symbol: share of the 50 single-letter
gene substitutions in its codon that give a working machine (synonymous or works). Then the
outcome of every single-letter deletion and insertion in the coded part, and of duplicating each
coded segment of 2, 4, 6, 8 letters (whole codons, in frame) or 1, 3 letters (out of frame)."""
import os, collections
g = {"__name__": "lib"}
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gm_scan.py")).read()
src = src.split('res = collections.Counter()')[0]
exec(src, g)
react, G, M, start, L = g["react"], g["G"], g["M"], g["start"], g["L"]

def outcome(g2):
    p = react(M, g2); ms = [x for x in p if x != g2]
    if not ms: return "no product", None
    if M in ms: return "synonymous", M
    m2 = ms[0]; p2 = react(m2, g2)
    return ("works" if (g2 in p2 and m2 in p2) else "broken"), m2

print("machine symbol: share of substitutions in its codon that keep a working machine")
row = []
for k in range(len(M)):
    ok = 0
    for j in (0, 1):
        i = start + 2 * k + j
        for c in L:
            if c == G[i]: continue
            o, _ = outcome(G[:i] + c + G[i + 1:]); ok += o in ("synonymous", "works")
    row.append(f"{M[k]}:{ok / 50:.2f}")
print("  " + "  ".join(row))
core = [i for i, ch in enumerate(M) if not ch.isalpha()]
print(f"  opcode positions {len(core)}, letter positions {len(M) - len(core)}")

def tally(lab, muts):
    t = collections.Counter(); ex = []
    for g2 in muts:
        o, m2 = outcome(g2); t[o] += 1
        if o == "works" and len(ex) < 3: ex.append(m2)
    n = sum(t.values())
    print(f"{lab:34} n={n:4}  " + "  ".join(f"{k} {t[k] / n:.2f}" for k in ("synonymous", "works", "broken", "no product")))
    for m2 in ex: print(f"      e.g. {m2}")
tally("single deletions", [G[:i] + G[i + 1:] for i in range(start, len(G))])
tally("single insertions", [G[:i] + c + G[i:] for i in range(start, len(G)) for c in "AEJOT"])
for w in (2, 4, 6, 8):
    tally(f"in-frame duplication of {w} letters", [G[:i + w] + G[i:i + w] + G[i + w:] for i in range(start, len(G) - w + 1, 2)])
for w in (1, 3):
    tally(f"out-of-frame duplication of {w}", [G[:i + w] + G[i:i + w] + G[i + w:] for i in range(start, len(G) - w + 1)])
