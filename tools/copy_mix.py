#!/usr/bin/env python3
"""Who copies whom over a long run: copy_mix.py RUN_DIR [--every 1000000] [--window 20000]

For the WINDOW steps before each multiple of EVERY, over exact copies (product == template),
using the placed column when present and otherwise attempts (older runs):
  self     maker and template the same species
  mutual   copies inside species pairs that copy each other (each direction >= 20)
  oneway   copies of a template by another species that the template does not copy back
  codeless share of copies whose template has no '=' (no copy loop of its own)
  pairs    effective number of mutual species pairs (inverse Simpson over pair totals)
"""
import argparse, collections, gzip, os

ap = argparse.ArgumentParser(); ap.add_argument("run")
ap.add_argument("--every", type=int, default=1000000); ap.add_argument("--window", type=int, default=20000)
a = ap.parse_args()
ropen = lambda p: open(p) if os.path.exists(p) else gzip.open(p + ".gz", "rt")
win = collections.defaultdict(collections.Counter); placed_col = None
for l in ropen(f"{a.run}/copies.tsv"):
    r = l.rstrip("\n").split("\t")
    t = int(r[0]); k = -(-t // a.every) * a.every
    if k - t >= a.window or r[2] != r[3]:
        continue
    if placed_col is None:
        placed_col = len(r) > 6
    win[k][(int(r[1]), int(r[2]))] += int(r[6]) if placed_col else int(r[4])
need = {t for c in win.values() for (m, t) in c}
seq = {}
for l in ropen(f"{a.run}/species.tsv"):
    i = l.find("\t")
    if i > 0 and l[:i].isdigit() and int(l[:i]) in need:
        seq[int(l[:i])] = l.rstrip("\n").split("\t")[-1]
print(f"== {a.run} ({'placed' if placed_col else 'attempted'} exact copies)")
print(f"{'step':>9} {'copies':>9} {'self':>5} {'mutual':>6} {'oneway':>6} {'codeless':>8} {'pairs':>6}")
for k in sorted(win):
    c = win[k]; tot = sum(c.values())
    if not tot:
        continue
    self_ = sum(v for (m, t), v in c.items() if m == t)
    pairs = {}
    for (m, t), v in c.items():
        if m < t and v >= 20 and c.get((t, m), 0) >= 20:
            pairs[(m, t)] = v + c[(t, m)]
    mut = sum(pairs.values())
    eff = (mut * mut / sum(v * v for v in pairs.values())) if mut else 0
    codeless = sum(v for (m, t), v in c.items() if "=" not in seq.get(t, "="))
    print(f"{k:>9} {tot:>9} {self_/tot:>5.2f} {mut/tot:>6.2f} {(tot-self_-mut)/tot:>6.2f} {codeless/tot:>8.2f} {eff:>6.1f}")
