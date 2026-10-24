#!/usr/bin/env python3
"""Summarise transplant scans: tp_agg.py SCAN_DIR

Files are <recipient>_<treatment>_rep<k>.json from transition_scan.py. Per recipient x
treatment: runs alive at the end, runs with processive share >= 0.10 at the end, median final
processive share, and detection / establishment steps (relative to the introduction).
"""
import collections, glob, json, os, statistics, sys

d = sys.argv[1]
g = collections.defaultdict(list)
for f in sorted(glob.glob(f"{d}/*.json")):
    b = os.path.basename(f)[:-5]
    if b.startswith("cache"):
        continue
    rec, tr, rep = b.rsplit("_", 2)
    o = json.load(open(f))
    rows = o["rows"]
    t0 = 0 if rec == "empty" else 2500000
    alive = bool(rows) and rows[-1]["n"] > 0 and rows[-1]["step"] >= t0 + 300000
    fin = rows[-1]["proc_share"] if rows else 0
    rel = lambda x: None if x is None else (x - t0) // 1000
    g[(rec, tr)].append((alive, fin, rel(o["detection"]), rel(o["establishment"])))
trs = ["E1", "E2", "M1", "M2", "ORD", "NONE"]
print("recipient x treatment: alive, established (>=10% processive at end), median final processive share,")
print("detection / establishment in k steps after introduction")
for rec in sorted({k[0] for k in g}, key=lambda r: (r != "empty", r)):
    print(f"== {rec}")
    for tr in trs:
        v = g.get((rec, tr))
        if not v:
            continue
        al = [x for x in v if x[0]]
        est = sum(1 for x in al if x[1] >= 0.10)
        med = statistics.median(x[1] for x in al) if al else 0
        print(f"  {tr:<5} alive {len(al)}/{len(v)}  established {est}/{len(v)}  median share {med:.2f}  det {[x[2] for x in v]}  est {[x[3] for x in v]}")
