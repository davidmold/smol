#!/usr/bin/env python3
"""gm_scan3.py -> can the machine's working core be duplicated? In-frame duplications (whole
codons) of 2-30 letters starting inside the core (machine symbols from the first opcode on),
outcome as gm_scan2.py. The ordinary-replicator route to self-scan is a duplication of the copy
loop; this asks whether the gene-and-machine code allows the same."""
import os, collections
g = {"__name__": "lib", "__file__": os.path.join(os.path.dirname(os.path.abspath(__file__)), "gm_scan2.py")}
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gm_scan2.py")).read()
src = src.split('print("machine symbol:')[0]
exec(src, g)
G, M, start, outcome = g["G"], g["M"], g["start"], g["outcome"]
first_op = next(i for i, ch in enumerate(M) if not ch.isalpha())
cs = start + 2 * first_op
print(f"core starts at machine symbol {first_op} ({M[first_op:]}), gene letter {cs}")
for w in (2, 4, 6, 8, 10, 12, 16, 20, 30):
    t = collections.Counter(); ex = []
    for i in range(cs, len(G) - w + 1, 2):
        o, m2 = outcome(G[:i + w] + G[i:i + w] + G[i + w:]); t[o] += 1
        if o == "works" and len(ex) < 2: ex.append(m2)
    n = sum(t.values())
    print(f"core in-frame duplication of {w:2} letters  n={n:3}  works {t['works'] / n:.2f}  broken {t['broken'] / n:.2f}  " + "  ".join(ex))
