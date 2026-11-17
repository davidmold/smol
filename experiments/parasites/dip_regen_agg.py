#!/usr/bin/env python3
"""dip_regen_agg.py DIR: parasite-free worlds (PARASITEFREE runs at 2M) continued 200k steps with
parasites allowed again (pf0m1), with parasites allowed but mutation off (pf0m0), and still
parasite-free (pf1m1). Population over time (median of runs), and for pf0m1 the step at which the
first parasite species (copied exactly >= 20 times, never a copier) appears."""
import sys, re, glob, collections, statistics
D = sys.argv[1]
offs = [0, 2000, 5000, 10000, 20000, 50000, 100000, 150000, 200000]
pop = collections.defaultdict(list); first = []
for f in glob.glob(f"{D}/w*/*/s*/log.txt"):
    c = f.split("/")[-3]
    vals = {int(m.group(1)): int(m.group(2)) for m in (re.match(r"Step (\d+) done, number of molecules is (\d+)", l) for l in open(f)) if m}
    base = 2000000
    vals.setdefault(base, max(vals.values()) if vals else 0)
    for o in offs:
        if base + o in vals: pop[(c, o)].append(vals[base + o])
    if c == "pf0m1":
        d = f.rsplit("/", 1)[0]
        copier, copied, when = set(), collections.Counter(), {}
        for l in open(f"{d}/copies.tsv"):
            p = l.split("\t")
            if len(p) < 5 or not p[0].isdigit(): continue
            st, a, pa, cc, n = int(p[0]), int(p[1]), int(p[2]), int(p[3]), int(p[4])
            if st <= base: continue
            if cc == pa or cc == a:
                copier.add(a); copied[cc] += n
                if copied[cc] >= 20 and cc not in when: when[cc] = st
        ps = [w for s, w in when.items() if s not in copier]
        first.append(min(ps) - base if ps else None)
print("population, median over 3 worlds x 5 seeds (n)")
for c in ("pf1m1", "pf0m0", "pf0m1"):
    print(f"  {c}: " + "  ".join(f"+{o//1000}k {statistics.median(pop[(c,o)]):.0f}" for o in offs if pop[(c, o)]))
f = [x for x in first if x is not None]
print(f"pf0m1: steps until the first successful parasite species: median {statistics.median(f):.0f}, range {min(f)}-{max(f)} ({len(f)}/{len(first)} runs)" if f else "no parasites")
