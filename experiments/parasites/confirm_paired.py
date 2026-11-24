#!/usr/bin/env python3
"""confirm_paired.py -> does escalation decline after 2M? Paired, per run: rho_run over 0.2-2M
(ts_gmt / ts_gmq) against rho_run over 2-6M (ts_gmt6 / ts_gmq6, the same runs continued).
Exploratory (suggested by an outside review after the results; not in CONFIRM_PLAN.md).
Wilcoxon signed-rank on the paired differences, median difference with a bootstrap 95% CI."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import confirm_stats as CS  # re-execs under the scipy venv if needed
import numpy as np, statistics as S
from scipy import stats
for early, late in (("gmt", "gmt6"), ("gmq", "gmq6")):
    a, _ = CS.analyse(early, CS.CC.check(early)); b, _ = CS.analyse(late, CS.CC.check(late))
    runs = sorted(set(a) & set(b), key=lambda x: int(x[1:]))
    e = np.array([a[r]["rho"] for r in runs]); l = np.array([b[r]["rho"] for r in runs]); d = l - e
    res = stats.bootstrap((d,), np.median, n_resamples=9999, method="percentile", random_state=CS.rng)
    print(f"{early}: {len(runs)} runs; median rho_run 0.2-2M {S.median(e):+.3f}, 2-6M {S.median(l):+.3f}; "
          f"late minus early median {np.median(d):+.3f} [95% CI {res.confidence_interval.low:+.3f}, "
          f"{res.confidence_interval.high:+.3f}]; less negative late in {int((d > 0).sum())}/{len(runs)}; "
          f"Wilcoxon p={stats.wilcoxon(d).pvalue:.2g}")
    for r, x, y in zip(runs, e, l): print(f"    {r}: {x:+.3f} -> {y:+.3f}")
