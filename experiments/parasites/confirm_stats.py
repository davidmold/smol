#!/usr/bin/env python3
"""confirm_stats.py SET [SET...] [--allow-incomplete] -> run-level tests of escalation on
time-shift matrices (SET or ts_SET, e.g. cf_gmt; matrices in ~/smol/runs/ts_SET).

Before scoring, confirm_check.py checks each set; a set with any technical gap is not scored
unless --allow-incomplete is given (then the gaps are listed with the results).

Data rules (CONFIRM_PLAN.md and its addendum):
- a cell counts only if both assay seeds are present (an assay extinction counts as population 0);
- runs that went extinct in evolution before the last matrix time are excluded and reported;
- d(H,P) = (pop(H,P) - pop(H,H)) / pop(H,H), mean of the two seeds, population at +10k steps.

Per run (the unit):
1. rho_run: Spearman correlation of d with the lag P - H over the run's off-diagonal cells
   (escalation: negative); rho_late on cells with H, P >= 0.6M (sets whose times start at 0.2M).
2. beta_H, beta_P: least squares of log(1 + pop(H,P)) on host time and parasite time (in M),
   all cells including the diagonal. beta_P < 0: later parasites hurt any fixed host more;
   beta_H > 0: later host worlds do better against any fixed parasites. Indirect evidence only:
   later worlds differ in more than their machines (see mech2.py for the direct assay).
Across runs: median with a bootstrap 95% CI, number of runs with the escalation sign, exact
binomial sign test, Wilcoxon signed-rank (scipy, exact where possible). Between the first set and
each other: Mann-Whitney U (scipy), the difference of medians with a bootstrap 95% CI, and the
rank-biserial correlation as effect size.
"""
import sys, os
try:
    import scipy  # noqa
except ImportError:  # re-run under the analysis venv, which has scipy
    py = os.path.expanduser("~/smenv/bin/python")
    if os.path.exists(py) and sys.executable != py: os.execv(py, [py] + sys.argv)
    raise
import re, glob, math, collections, statistics as S
import numpy as np
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import confirm_check as CC
R = CC.R; OFF = CC.OFF
rng = np.random.default_rng(1)

def load(name):
    seeds = collections.defaultdict(dict)
    for f in glob.glob(f"{R}/ts_{name}/*/h*_p*/s*/log.txt"):
        parts = f.split("/"); r, cell, sd = parts[-4], parts[-3], parts[-2]
        m = re.match(r"h(\d+)k_p(\d+)k", cell); h, p = int(m.group(1)) * 1000, int(m.group(2)) * 1000
        o = CC.swap_outcome(f, h)
        if o == "assay_extinct": v = 0
        elif o == "ok":
            v = int(re.search(rf"^Step {h + OFF} done, number of molecules is (\d+)", open(f, errors="replace").read(), re.M).group(1))
        else: continue
        seeds[(r, h, p)][sd] = v
    return {k: sum(v.values()) / 2 for k, v in seeds.items() if len(v) == 2}

def ols2(H, P, Y):
    X = np.column_stack([np.ones(len(H)), H, P]); b = np.linalg.lstsq(X, np.array(Y), rcond=None)[0]
    return b[1], b[2]

def analyse(name, chk):
    mean = load(name); times = CC.EXPECT[name][1]
    excl = {r for r, v in chk["runs"].items() if v["extinct"] is not None and v["extinct"] < times[-1]}
    out = {}
    for r in sorted({k[0] for k in mean} - excl, key=lambda x: int(x[1:])):
        cells = {(h, p): v for (rr, h, p), v in mean.items() if rr == r}
        d = {(h, p): (v - cells[(h, h)]) / cells[(h, h)] for (h, p), v in cells.items() if p != h and cells.get((h, h))}
        o = {"cells": len(d)}
        if len(d) >= 3:
            o["rho"] = stats.spearmanr([p - h for h, p in d], list(d.values())).statistic
            if times[0] < 600000:
                late = {k: x for k, x in d.items() if k[0] >= 600000 and k[1] >= 600000}
                if len(late) >= 3: o["rho_late"] = stats.spearmanr([p - h for h, p in late], list(late.values())).statistic
        if len(cells) >= 4:
            o["bH"], o["bP"] = ols2([h / 1e6 for h, p in cells], [p / 1e6 for h, p in cells], [math.log1p(v) for v in cells.values()])
        out[r] = o
    return out, excl

def ci_median(v):
    if len(v) < 3: return (float("nan"), float("nan"))
    res = stats.bootstrap((np.array(v),), np.median, n_resamples=9999, method="percentile", random_state=rng)
    return res.confidence_interval.low, res.confidence_interval.high

def report(name, out, excl, chk):
    print(f"\n== {name}: {len(out)} runs scored; excluded (extinct in evolution): {sorted(excl) or 'none'}")
    print(f"   off-diagonal cells per run: min {min((o['cells'] for o in out.values()), default=0)}, "
          f"max {max((o['cells'] for o in out.values()), default=0)}; windows without parasites: "
          f"{chk['count'].get('no_parasite_windows', 0)}; assay extinctions: {chk['count'].get('swaps_assay_extinct', 0)}")
    for key, lab, esc in (("rho", "rho(lag, d), all cells", -1), ("rho_late", "rho(lag, d), H,P >= 0.6M", -1),
                          ("bP", "beta_P (parasite time)", -1), ("bH", "beta_H (host time)", +1)):
        v = [o[key] for o in out.values() if key in o and not math.isnan(o[key])]
        if len(v) < 3: continue
        k = sum(1 for x in v if x * esc > 0); lo, hi = ci_median(v)
        nz = [x for x in v if x != 0]
        wp = stats.wilcoxon(nz).pvalue if len(nz) >= 3 else float("nan")
        print(f"  {lab:28} median {S.median(v):+.3f} [95% CI {lo:+.3f}, {hi:+.3f}]  escalation sign {k}/{len(v)}"
              f"  sign p={stats.binomtest(k, len(v)).pvalue:.2g}  Wilcoxon p={wp:.2g}")
    both = sum(1 for o in out.values() if o.get("bP", 0) < 0 and o.get("bH", 0) > 0)
    print(f"  runs with beta_P < 0 and beta_H > 0: {both}/{len(out)}")
    for r, o in out.items():
        print(f"    {r}: cells {o['cells']}  " + "  ".join(f"{k} {o[k]:+.3f}" for k in ("rho", "rho_late", "bP", "bH") if k in o))

def compare(n1, a, n2, b):
    if len(a) < 3 or len(b) < 3: return
    u = stats.mannwhitneyu(a, b, alternative="two-sided")
    rb = 1 - 2 * u.statistic / (len(a) * len(b))       # rank-biserial: > 0 means a tends lower
    diff = lambda x, y: np.median(x) - np.median(y)
    res = stats.bootstrap((np.array(a), np.array(b)), diff, n_resamples=9999, method="percentile", random_state=rng)
    print(f"\nrho_run, {n1} vs {n2}: median difference {S.median(a) - S.median(b):+.3f} "
          f"[95% CI {res.confidence_interval.low:+.3f}, {res.confidence_interval.high:+.3f}]  "
          f"Mann-Whitney p={u.pvalue:.2g}  rank-biserial {rb:+.2f} (> 0: {n1} lower)")

if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    allow = "--allow-incomplete" in sys.argv
    names = [a[3:] if a.startswith("ts_") else a for a in args]
    res = []; failed = False
    print("Completeness check (confirm_check.py):")
    for n in names:
        chk = CC.check(n); CC.report(n, chk)
        if chk["technical"] and not allow:
            failed = True; continue
        out, excl = analyse(n, chk); res.append((n, out, excl, chk))
    if failed:
        print("\nNOT SCORED: technical gaps above. Fix and rerun the missing jobs, or score with "
              "--allow-incomplete and report the gaps.")
    for n, out, excl, chk in res: report(n, out, excl, chk)
    if len(res) > 1:
        a = [o["rho"] for o in res[0][1].values() if "rho" in o]
        for n, out, _, _ in res[1:]:
            compare(res[0][0], a, n, [o["rho"] for o in out.values() if "rho" in o])
    sys.exit(2 if failed else 0)
