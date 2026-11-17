#!/usr/bin/env python3
"""dropone.py setup OUT / score OUT: run 4's two puzzling pools (src60000, from before D existed,
and src100000, the original competition pool) rerun with each parasite left out in turn (the rest
rescaled to 1000), 10 seeds each, same set-up as ratio.py. Which parasite is responsible for D losing?"""
import sys, os, glob, random, statistics
H = os.path.expanduser("~"); mode, OUT = sys.argv[1], os.path.abspath(sys.argv[2])
import importlib.util
spec = importlib.util.spec_from_file_location("ko", os.path.join(os.path.dirname(os.path.abspath(__file__)), "knockout.py"))
sys.argv = [sys.argv[0], "lib", OUT]; ko = importlib.util.module_from_spec(spec); spec.loader.exec_module(ko)
POOLS = [("run4", "src60000"), ("run4", "src100000")]
if mode == "setup":
    jobs = open(f"{OUT}/jobs.txt", "w")
    for run, nm in POOLS:
        par, cnt, _ = ko.P[(run, nm)]; A, D = ko.PAIRS[run][2], ko.PAIRS[run][3]
        for drop in range(len(par)):
            keep = [(s, n) for i, (s, n) in enumerate(zip(par, cnt)) if i != drop]; t = sum(n for _, n in keep)
            for cs in range(1, 11):
                d = f"{OUT}/{run}/{nm}/drop{drop}/s{cs}"; os.makedirs(d, exist_ok=True)
                rng = random.Random(cs * 1000 + int(run[3:]))
                mix = [(A, 1000), (D, 1000)] + [(s, round(1000 * n / t)) for s, n in keep]
                cells = rng.sample([(x, y) for x in range(125) for y in range(100)], sum(n for _, n in mix))
                lines = [f"RANDSEED {cs}", "NSTEPS 30000", "MUTATE 0", "COMPAT_DELETION 1",
                         "SNAPEVERY 100000000", "CKPTEVERY 1000000", "POPDYEVERY 500"] + ko.BASE
                k = 0
                for s, n in mix:
                    lines.append(f"AGENT {s} {n} Q")
                    for _ in range(n):
                        x, y = cells[k]; k += 1; lines.append(f"GRIDPOS {x} {y}")
                open(f"{d}/run.conf", "w").write("\n".join(lines) + "\n")
                jobs.write(f"cd {d} && {ko.SMOL} run run.conf > log.txt 2>&1\n")
else:
    for run, nm in POOLS:
        par, cnt, full = ko.P[(run, nm)]; A, D = ko.PAIRS[run][2], ko.PAIRS[run][3]
        print(f"== {run} {nm}: full pool D share {ko.med(glob.glob(full + '/s*'), A, D):.2f}")
        for drop, (s, n) in enumerate(zip(par, cnt)):
            print(f"  without {s:<62} ({n:>4} in source)  D share {ko.med(glob.glob(f'{OUT}/{run}/{nm}/drop{drop}/s*'), A, D):.2f}")
