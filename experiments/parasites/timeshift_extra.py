#!/usr/bin/env python3
"""timeshift_extra.py DIR: extra tests on a time-shift matrix (timeshift_agg.py's d values):
rank correlation of d with the lag (escalation: negative), for all hosts and for subsets, and the
median d of the parasites from each time against all other host times (negative = those
parasites are more harmful than each host's own)."""
import sys, math, statistics as S
D = sys.argv[1]
src = open(__file__.replace("timeshift_extra.py", "timeshift_agg.py")).read().split("step = times[1]")[0]
src = src.replace("D = sys.argv[1]", "D = %r" % D).replace("int(sys.argv[2]) if len(sys.argv) > 2 else 10000", "10000")
g = {"__file__": __file__}; exec(src, g); d = g["d"]
def sp(x, y):
    def rk(v):
        s = sorted(v); return [(s.index(a) + len(s) - 1 - s[::-1].index(a)) / 2 for a in v]
    a, b = rk(x), rk(y); n = len(x); ma, mb = sum(a) / n, sum(b) / n
    return sum((i - ma) * (j - mb) for i, j in zip(a, b)) / math.sqrt(sum((i - ma) ** 2 for i in a) * sum((j - mb) ** 2 for j in b))
cells = [(r, h, p, x) for (r, h, p), x in d.items() if p != h]
def summ(lab, cs):
    rho = sp([p - h for r, h, p, x in cs], [x for *_, x in cs]); n = len(cs); t = rho * math.sqrt((n - 2) / (1 - rho * rho))
    fut = [x for r, h, p, x in cs if p > h]; past = [x for r, h, p, x in cs if p < h]
    print(f"{lab:44} n={n:4}  rho(lag, d) = {rho:+.2f} (t = {t:+.1f})  past {S.median(past):+.3f}  future {S.median(fut):+.3f}")
summ("all hosts", cells)
summ("hosts 0.6-1.8M", [c for c in cells if 600000 <= c[1] <= 1800000])
summ("hosts and parasites both from 0.6M on", [c for c in cells if c[1] >= 600000 and c[2] >= 600000])
print("\nParasites from time P against all other host times: median d (< 0: more harmful than the host's own)")
for p in sorted({c[2] for c in cells}):
    v = [x for r, h, pp, x in cells if pp == p]
    print(f"  P = {p / 1e6:.1f}M   {S.median(v):+.3f}   (< 0 in {sum(1 for x in v if x < 0)}/{len(v)})")
