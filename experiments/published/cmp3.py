#!/usr/bin/env python3
"""cmp3.py: survivors of the paper's runs, our Oct rerun (old rules) and the COMPAT_DELETION rerun,
compared on the assay (tools/assay.py) averaged over 1.0-1.2M and 1.9-2.0M. Mann-Whitney U
(two-sided, normal approximation with tie correction) between the paper and each rerun."""
import sys, statistics, math, collections
def load(fn):
    runs = {}; cur = None
    for l in open(fn):
        if l.startswith("=="): cur = l.split()[1]; runs[cur] = {}; continue
        p = l.split()
        if len(p) == 9 and p[0].isdigit(): runs[cur][int(p[0])] = [float(x) for x in p]
    return {k: v for k, v in runs.items() if 2000000 in v}
def win(r, lo, hi, k):
    v = [r[s][k] for s in r if lo <= s <= hi and not math.isnan(r[s][k])]
    return statistics.mean(v) if v else float("nan")
def mwu(a, b):
    allv = sorted([(x, 0) for x in a] + [(x, 1) for x in b]); n1, n2 = len(a), len(b)
    ranks = {}; i = 0
    while i < len(allv):
        j = i
        while j < len(allv) and allv[j][0] == allv[i][0]: j += 1
        for k in range(i, j): ranks.setdefault(k, (i + j + 1) / 2)
        i = j
    R1 = sum(ranks[k] for k, (x, g) in enumerate(allv) if g == 0)
    U = R1 - n1 * (n1 + 1) / 2; mu = n1 * n2 / 2; sd = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    z = (U - mu) / sd if sd else 0
    return 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
sets = {"paper (8/20)": load("assay_pub.txt"), "Oct rerun, old rules (5/15)": load("assay_pubrerun.txt"),
        "rerun + COMPAT_DELETION (7/24)": load("assay_pubdel.txt")}
cols = {"N": 1, "selfrep": 3, "mutual": 4, "parasitic": 5, "rtime": 6, "toplen": 7}
for lo, hi in ((1000000, 1200000), (1900000, 2000000)):
    print(f"== mean over {lo // 1000}k-{hi // 1000}k, median across surviving runs [range]; p = Mann-Whitney vs paper")
    for name, runs in sets.items():
        row = []
        for c, k in cols.items():
            v = [win(r, lo, hi, k) for r in runs.values()]
            row.append(f"{c} {statistics.median(v):7.2f} [{min(v):.2f}-{max(v):.2f}]" if c not in ("N", "rtime", "toplen")
                       else f"{c} {statistics.median(v):6.0f} [{min(v):.0f}-{max(v):.0f}]")
        print(f"  {name:32} n={len(runs)}  " + "  ".join(row))
    base = sets["paper (8/20)"]
    for name, runs in list(sets.items())[1:]:
        ps = []
        for c, k in cols.items():
            ps.append(f"{c} p={mwu([win(r, lo, hi, k) for r in base.values()], [win(r, lo, hi, k) for r in runs.values()]):.2f}")
        print(f"  {name:32} vs paper: " + "  ".join(ps))
# change from 1.0-1.2M to 1.9-2.0M within each run
print("== change from 1.0-1.2M to 1.9-2.0M within each surviving run (median [range])")
for name, runs in sets.items():
    row = []
    for c in ("parasitic", "mutual", "rtime", "N"):
        k = cols[c]; d = [win(r, 1900000, 2000000, k) - win(r, 1000000, 1200000, k) for r in runs.values()]
        row.append(f"{c} {statistics.median(d):+.2f} [{min(d):+.2f},{max(d):+.2f}]" if c in ("parasitic", "mutual") else f"{c} {statistics.median(d):+.0f}")
    print(f"  {name:32} " + "  ".join(row))
