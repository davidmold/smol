#!/usr/bin/env python3
"""Score the disturbance experiment (disturb.py; run on aida).

    disturb_agg.py [--scan] [--end STEP]

--scan: run transition_scan.py --sample 1000 --every 10000 on every finished run without a
scan.json (shared cache). --end: last step to score (default 9.9M + 200k; the continued
populations are scored to 400k separately with --end 10300000 and --tag cont).
Per run: any processive species seen in a sample; detection (a sampled processive species with
>= 20 molecules); established (share >= 10% for >= 3 consecutive samples, the scan's rule); share
at 9.9M+150k (disturbance stops) and at the end; extinct; mean molecule count; placed copies;
molecules removed by disturbance. Per treatment: counts, medians, and Fisher exact tests on
established runs between treatments. Established runs' top processive species are listed.
"""
import glob, json, os, re, statistics as st, subprocess, sys
from math import comb

H = os.path.expanduser("~/smol"); R = f"{H}/runs/disturb"; CACHE = "/tmp/c_disturb.json"
T0 = 9900000; STOP = T0 + 150000
end = int(sys.argv[sys.argv.index("--end") + 1]) if "--end" in sys.argv else T0 + 200000
tag = sys.argv[sys.argv.index("--tag") + 1] if "--tag" in sys.argv else ""
SJ = f"scan{tag}.json"

def done(d):
    t = open(f"{d}/log.txt", errors="replace").read()[-3000:]
    return f"Step {end} done" in t or "EXTINCT" in t

runs = sorted(glob.glob(f"{R}/*/rep*"), key=lambda d: (d.split("/")[-2], int(d.split("rep")[-1])))
if "--scan" in sys.argv:
    if not os.path.exists(CACHE):
        for c in ("/tmp/c_rep.json", "/tmp/c_natural.json", "/tmp/cs2.json"):
            if os.path.exists(c):
                subprocess.run(["cp", c, CACHE]); break
    for d in runs:
        if done(d) and not os.path.exists(f"{d}/{SJ}"):
            subprocess.run([sys.executable, "-I", f"{H}/tools/transition_scan.py", d, "--start", str(T0),
                            "--max", str(end), "--every", "10000", "--sample", "1000", "--cache", CACHE,
                            "--json", f"{d}/{SJ}"], stdout=open(f"{d}/scan{tag}.txt", "w"), check=True)

def fisher_two_sided(a, n1, b, n2):
    """Two-sided Fisher exact p for a/n1 vs b/n2."""
    k = a + b; N = n1 + n2
    pr = lambda x: comb(n1, x) * comb(n2, k - x) / comb(N, k)
    p0 = pr(a)
    return min(1.0, sum(pr(x) for x in range(max(0, k - n2), min(k, n1) + 1) if pr(x) <= p0 * (1 + 1e-9)))

res = {}
for d in runs:
    if not os.path.exists(f"{d}/{SJ}"):
        continue
    t = d.split("/")[-2]
    J = json.load(open(f"{d}/{SJ}")); rows = [r for r in J["rows"] if r["step"] <= end]
    log = open(f"{d}/log.txt", errors="replace").read()
    run_len = 0; est = None
    for i, r in enumerate(rows):
        run_len = run_len + 1 if r["proc_share"] >= 0.10 else 0
        if est is None and run_len >= 3:
            est = rows[i - 2]["step"]
    at = lambda s: next((r["proc_share"] for r in rows if r["step"] == s), None)
    res.setdefault(t, []).append(dict(
        run=d, seen=any(r["proc_species"] > 0 for r in rows), det=J["detection"], est=est,
        stop=at(STOP), end=rows[-1]["proc_share"] if rows else 0, ext="EXTINCT" in log,
        n=st.mean(r["n"] for r in rows) if rows else 0, cum=rows[-1]["cum"] if rows else 0,
        removed=sum(int(x) for x in re.findall(r"DISTURB to step \d+: removed (\d+)", log))))

print(f"scored to step {end}")
print(f"{'treatment':10} {'runs':>4} {'seen':>4} {'det':>4} {'estab':>5} {'end>=10%':>8} {'kept@end':>8} {'extinct':>7} "
      f"{'mean_N':>7} {'placed_copies':>13} {'removed':>8}")
for t in ("base", "clustered", "dispersed"):
    v = res.get(t, [])
    if not v:
        continue
    print(f"{t:10} {len(v):>4} {sum(x['seen'] for x in v):>4} {sum(x['det'] is not None for x in v):>4} "
          f"{sum(x['est'] is not None for x in v):>5} {sum(x['end'] >= 0.10 for x in v):>8} "
          f"{sum(x['est'] is not None and x['end'] >= 0.10 for x in v):>8} {sum(x['ext'] for x in v):>7} "
          f"{st.median(x['n'] for x in v):>7.0f} {st.median(x['cum'] for x in v):>13.3g} {st.median(x['removed'] for x in v):>8.0f}")
print()
E = {t: (sum(x["est"] is not None for x in v), len(v)) for t, v in res.items()}
for a_, b_ in (("clustered", "base"), ("dispersed", "base"), ("clustered", "dispersed")):
    if a_ in E and b_ in E:
        print(f"established {a_} {E[a_][0]}/{E[a_][1]} vs {b_} {E[b_][0]}/{E[b_][1]}: "
              f"two-sided Fisher p = {fisher_two_sided(E[a_][0], E[a_][1], E[b_][0], E[b_][1]):.3f}")
print()
for t, v in res.items():
    for x in v:
        if x["est"] is not None or x["end"] >= 0.10:
            print(f"{t:10} {x['run'].split('/')[-1]:6} est {x['est']} det {x['det']} share@stop {x['stop']} "
                  f"share@end {x['end']:.3f}")
