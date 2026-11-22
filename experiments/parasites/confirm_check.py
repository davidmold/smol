#!/usr/bin/env python3
"""confirm_check.py SET [SET...] -> completeness check of evolution runs and time-shift matrices,
run before any scoring (confirm_stats.py calls it and refuses to score a set that fails).

For each set (run dirs in ~/smol/runs/SET/sN, matrix in ~/smol/runs/ts_SET) and each expected
run, time and cell, every gap is put in one of these classes:
  evo_extinct     the evolution run went extinct (its log says EXTINCT) before the last time;
                  times after it are expected to be missing
  no_parasites    the classification window at that time has no gene-parasites, so the
                  time-shift script makes no cells with it as host or parasite time (by design)
  assay_extinct   a swap run went extinct by +OFF steps: a valid outcome (population 0), counted
  technical       anything else missing: a run that stopped early without EXTINCT, a missing
                  checkpoint, window, classes.tsv, swap directory or log, or a log without the
                  +OFF line. These make the set fail.
Prints a table per set and, with --json F, writes the details. Exit status 1 if any set has a
technical gap.
"""
import sys, os, re, json, glob, collections
R = os.path.expanduser("~/smol/runs"); OFF = 10000; SEEDS = (1, 2)
T2 = [200000 * i for i in range(1, 11)]
T6 = [2000000 + 500000 * i for i in range(9)]
EXPECT = {"gmt": (range(301, 311), T2), "gmq": (range(301, 311), T2), "gm": (range(301, 311), T2),
          "cf_gmt": (range(311, 341), T2), "cf_gmq": (range(311, 341), T2), "cf_gmd": (range(311, 341), T2),
          "gmt6": (range(301, 311), T6), "gmq6": (range(301, 311), T6)}

def run_status(d):
    """(last step reached, extinct step or None) from an evolution run's log."""
    last, ext = 0, None
    try:
        for l in open(f"{d}/log.txt", errors="replace"):
            m = re.match(r"Step (\d+) done", l)
            if m: last = int(m.group(1))
            if l.startswith("EXTINCT"):
                m = re.search(r"step (\d+)", l); ext = int(m.group(1)) if m else last
    except FileNotFoundError:
        return None, None
    return last, ext

def swap_outcome(f, h):
    try:
        txt = open(f, errors="replace").read()
    except FileNotFoundError:
        return "technical:no_log"
    if re.search(rf"^Step {h + OFF} done", txt, re.M): return "ok"
    m = re.search(r"^EXTINCT.*?step (\d+)", txt, re.M)
    if m and int(m.group(1)) <= h + OFF: return "assay_extinct"
    return "technical:no_line"

def check(name):
    seeds, times = EXPECT[name]
    S, TS = f"{R}/{name}", f"{R}/ts_{name}"
    rep = {"runs": {}, "count": collections.Counter(), "technical": []}
    for s in seeds:
        r = f"s{s}"; d = f"{S}/{r}"; last, ext = run_status(d)
        info = {"last": last, "extinct": ext, "times": {}}
        rep["runs"][r] = info
        if last is None:
            rep["technical"].append(f"{r}: no run log"); continue
        if ext is None and last < times[-1]:
            rep["technical"].append(f"{r}: stopped at {last} without EXTINCT")
        if ext is not None and ext < times[-1]:
            rep["count"]["evo_extinct_runs"] += 1
        alive = [t for t in times if ext is None or t < ext]
        withP = []
        for t in alive:
            if t > last: continue
            if not os.path.exists(f"{d}/checkpoint_{t}.txt.gz"):
                rep["technical"].append(f"{r}: no checkpoint_{t}"); continue
            cl = f"{TS}/{r}/w{t // 1000}k/classes.tsv"
            if not os.path.exists(cl) or os.path.getsize(cl) == 0:
                rep["technical"].append(f"{r}: no classes.tsv at {t}"); continue
            hasP = any(l.split()[1] == "P" for l in open(cl) if len(l.split()) >= 2)
            info["times"][t] = "P" if hasP else "noP"
            if hasP: withP.append(t)
            else: rep["count"]["no_parasite_windows"] += 1
        for h in withP:
            for p in withP:
                for sd in SEEDS:
                    o = swap_outcome(f"{TS}/{r}/h{h // 1000}k_p{p // 1000}k/s{sd}/log.txt", h)
                    if o.startswith("technical"):
                        rep["technical"].append(f"{r}: cell h{h // 1000}k_p{p // 1000}k seed {sd}: {o.split(':')[1]}")
                    else:
                        rep["count"]["swaps_" + o] += 1
    return rep

def report(name, rep):
    c = rep["count"]; n = len(rep["runs"])
    print(f"{name}: {n} runs expected; evo-extinct before the last time: {c['evo_extinct_runs']}; "
          f"windows without parasites: {c['no_parasite_windows']}; swaps ok: {c['swaps_ok']}, "
          f"assay-extinct: {c['swaps_assay_extinct']}; technical gaps: {len(rep['technical'])}")
    for t in rep["technical"][:20]: print("   TECHNICAL", t)
    if len(rep["technical"]) > 20: print(f"   ... and {len(rep['technical']) - 20} more")

if __name__ == "__main__":
    args = sys.argv[1:]; jf = None
    if "--json" in args:
        i = args.index("--json"); jf = args[i + 1]; del args[i:i + 2]
    out = {}; bad = False
    for name in args:
        rep = check(name); report(name, rep); bad |= bool(rep["technical"])
        out[name] = {"count": dict(rep["count"]), "technical": rep["technical"],
                     "runs": {r: {"last": v["last"], "extinct": v["extinct"],
                                  "times": {str(k): x for k, x in v["times"].items()}} for r, v in rep["runs"].items()}}
    if jf: json.dump(out, open(jf, "w"), indent=1)
    sys.exit(1 if bad else 0)
