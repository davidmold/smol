#!/usr/bin/env python3
"""timeshift_agg.py DIR: summarise the time-shift matrix (timeshift.sh).

For run r, host time H and parasite time P: pop(H,P) = population OFF steps after H's parasites
were replaced by parasites from time P (mean of 3 seeds). Effect of foreign parasites:
    d(H,P) = (pop(H,P) - pop(H,H)) / pop(H,H)
d > 0: the host world does better against parasites from P than against its own (P's parasites
are less harmful); d < 0: they are more harmful.

Readings (host-parasite coevolution, e.g. Gandon et al. 2008):
  escalation (arms race): past parasites less harmful (d > 0 for P < H), future parasites more
                          harmful (d < 0 for P > H): d falls with the lag P - H.
  fluctuating / Red Queen: own-time parasites most harmful; d > 0 on both sides of the diagonal,
                          possibly falling back at long lags (old parasites work again).
Prints: d by lag (median over runs and hosts, with the share of cells > 0), the median matrix,
and per run the past-vs-future comparison.
"""
import sys, os, re, glob, statistics as S, collections
D = sys.argv[1]
OFF = int(sys.argv[2]) if len(sys.argv) > 2 else 10000
pop = collections.defaultdict(list)
for f in glob.glob(f"{D}/*/h*_p*/s*/log.txt"):
    r, cell = f.split("/")[-4], f.split("/")[-3]
    m = re.match(r"h(\d+)k_p(\d+)k", cell); h, p = int(m.group(1)) * 1000, int(m.group(2)) * 1000
    v = None
    for l in open(f, errors="replace"):
        mm = re.match(r"Step (\d+) done, number of molecules is (\d+)", l)
        if mm and int(mm.group(1)) == h + OFF: v = int(mm.group(2))
        if l.startswith("EXTINCT"):
            st = int(re.search(r"step (\d+)", l).group(1))
            if st <= h + OFF: v = 0
    if v is not None: pop[(r, h, p)].append(v)
mean = {k: sum(v) / len(v) for k, v in pop.items() if v}
runs = sorted({k[0] for k in mean}, key=lambda x: int(x[1:]))
times = sorted({k[1] for k in mean})
d = {}
for (r, h, p), v in mean.items():
    base = mean.get((r, h, h))
    if base: d[(r, h, p)] = (v - base) / base
step = times[1] - times[0]
print(f"Relative population change at +{OFF} steps when a world's parasites come from another time\n"
      f"(d > 0: foreign parasites less harmful than its own). {len(runs)} runs, {len(times)} times, {len(mean)} cells.\n")
print("lag (parasite time - host time)   median d   cells with d > 0")
bylag = collections.defaultdict(list)
for (r, h, p), x in d.items():
    if p != h: bylag[(p - h) // step].append(x)
for lag in sorted(bylag):
    v = bylag[lag]
    print(f"  {lag * step / 1e6:+5.1f}M                          {S.median(v):+6.3f}     {sum(1 for x in v if x > 0)}/{len(v)}")
print("\nMedian d over runs (rows: host time, columns: parasite time, in M)")
print("host\\par " + "".join(f"{t / 1e6:7.1f}" for t in times))
for h in times:
    row = []
    for p in times:
        v = [d[(r, h, p)] for r in runs if (r, h, p) in d]
        row.append(f"{S.median(v):+7.2f}" if v and p != h else "      .")
    print(f"  {h / 1e6:4.1f}   " + "".join(row))
print("\nPer run: mean d for parasites from the past vs from the future (all hosts)")
for r in runs:
    past = [x for (rr, h, p), x in d.items() if rr == r and p < h]
    fut = [x for (rr, h, p), x in d.items() if rr == r and p > h]
    if past and fut:
        print(f"  {r}: past {S.mean(past):+.3f}  future {S.mean(fut):+.3f}")
past = [x for (r, h, p), x in d.items() if p < h]; fut = [x for (r, h, p), x in d.items() if p > h]
print(f"\nAll: past parasites d median {S.median(past):+.3f} (> 0 in {sum(1 for x in past if x > 0)}/{len(past)}); "
      f"future parasites d median {S.median(fut):+.3f} (> 0 in {sum(1 for x in fut if x > 0)}/{len(fut)})")
