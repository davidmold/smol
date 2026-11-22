#!/usr/bin/env python3
"""confirm_stats.py TS_DIR... -> run-level tests of escalation on time-shift matrices.

The pooled rank correlation in timeshift_extra.py treats every cell as independent, but the cells
of one run share a world (and each row shares its baseline pop(H,H)), so its t overstates the
evidence. Here each run gives one number and runs are the units:

1. rho_run: Spearman rank correlation of d(H,P) with the lag P - H over that run's off-diagonal
   cells (escalation: negative). Also on the late block (H, P >= 0.6M).
2. Two-way decomposition of raw outcomes: least squares of log(1 + pop(H,P)) on H and P (in M),
   all cells including the diagonal. An arms race needs both sides to move:
     beta_P < 0  later parasites hurt any fixed host more   (parasites improve)
     beta_H > 0  later hosts do better against any fixed parasites (hosts improve)
   Parasites improving against static hosts gives beta_P < 0 with beta_H ~ 0; hosts winning a
   one-sided race gives beta_H > 0 with beta_P ~ 0.
Across runs: median, count of runs with the escalation sign, and a two-sided sign test and
Wilcoxon signed-rank p (exact for n <= 30 via enumeration-free normal approx otherwise).
With two or more TS_DIRs, also Mann-Whitney between the first and each other (rho_run).
"""
import sys, os, re, glob, math, statistics as S, collections

def load(D, OFF=10000):
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
    return {k: sum(v) / len(v) for k, v in pop.items() if v}

def rank(v):
    s = sorted(v); return [(s.index(a) + len(s) - 1 - s[::-1].index(a)) / 2 for a in v]
def spear(x, y):
    a, b = rank(x), rank(y); n = len(x); ma, mb = sum(a) / n, sum(b) / n
    den = math.sqrt(sum((i - ma) ** 2 for i in a) * sum((j - mb) ** 2 for j in b))
    return sum((i - ma) * (j - mb) for i, j in zip(a, b)) / den if den else 0.0
def ols2(H, P, Y):
    n = len(Y); mh, mp, my = sum(H) / n, sum(P) / n, sum(Y) / n
    h = [x - mh for x in H]; p = [x - mp for x in P]; y = [x - my for x in Y]
    shh = sum(a * a for a in h); spp = sum(a * a for a in p); shp = sum(a * b for a, b in zip(h, p))
    shy = sum(a * b for a, b in zip(h, y)); spy = sum(a * b for a, b in zip(p, y))
    det = shh * spp - shp * shp
    return (spp * shy - shp * spy) / det, (shh * spy - shp * shy) / det
def phi(z): return 0.5 * (1 + math.erf(z / math.sqrt(2)))
def sign_p(k, n):
    m = min(k, n - k); return min(1.0, 2 * sum(math.comb(n, i) for i in range(m + 1)) / 2 ** n)
def wilcoxon_p(v):
    v = [x for x in v if x != 0]; n = len(v)
    if n < 5: return float("nan")
    r = rank([abs(x) for x in v]); w = sum(ri + 1 for ri, x in zip(r, v) if x > 0)
    mu = n * (n + 1) / 4; sd = math.sqrt(n * (n + 1) * (2 * n + 1) / 24)
    return 2 * (1 - phi(abs(w - mu) / sd))
def mannwhitney_p(a, b):
    r = rank(a + b); ra = sum(x + 1 for x in r[:len(a)]); n1, n2 = len(a), len(b)
    u = ra - n1 * (n1 + 1) / 2; mu = n1 * n2 / 2; sd = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    return 2 * (1 - phi(abs(u - mu) / sd))

def analyse(D):
    mean = load(D); runs = sorted({k[0] for k in mean}, key=lambda x: int(x[1:]))
    out = collections.defaultdict(dict)
    for r in runs:
        cells = {(h, p): v for (rr, h, p), v in mean.items() if rr == r}
        d = {(h, p): (v - cells[(h, h)]) / cells[(h, h)] for (h, p), v in cells.items()
             if p != h and cells.get((h, h))}
        if len(d) >= 10:
            out[r]["rho"] = spear([p - h for h, p in d], list(d.values()))
            late = {k: x for k, x in d.items() if k[0] >= 600000 and k[1] >= 600000}
            if len(late) >= 10: out[r]["rho_late"] = spear([p - h for h, p in late], list(late.values()))
        if len(cells) >= 10:
            bh, bp = ols2([h / 1e6 for h, p in cells], [p / 1e6 for h, p in cells],
                          [math.log1p(v) for v in cells.values()])
            out[r]["bH"], out[r]["bP"] = bh, bp
    return out

def report(D, out):
    print(f"\n== {D}: {len(out)} runs")
    for key, lab, esc in (("rho", "rho(lag, d), all cells", -1), ("rho_late", "rho(lag, d), H,P >= 0.6M", -1),
                          ("bP", "beta_P (parasite time)", -1), ("bH", "beta_H (host time)", +1)):
        v = [o[key] for o in out.values() if key in o]
        if not v: continue
        k = sum(1 for x in v if x * esc > 0)
        print(f"  {lab:28} median {S.median(v):+.3f}  escalation sign in {k}/{len(v)} runs"
              f"  sign p={sign_p(k, len(v)):.2g}  Wilcoxon p={wilcoxon_p(v):.2g}")
    both = sum(1 for o in out.values() if o.get("bP", 0) < 0 and o.get("bH", 0) > 0)
    print(f"  runs with both sides moving (beta_P < 0 and beta_H > 0): {both}/{len(out)}")
    for r, o in out.items():
        print(f"    {r}: " + "  ".join(f"{k} {o[k]:+.3f}" for k in ("rho", "rho_late", "bP", "bH") if k in o))

if __name__ == "__main__":
    res = [(D, analyse(D)) for D in sys.argv[1:]]
    for D, o in res: report(D, o)
    if len(res) > 1:
        a = [o["rho"] for o in res[0][1].values() if "rho" in o]
        for D, o in res[1:]:
            b = [x["rho"] for x in o.values() if "rho" in x]
            print(f"\nMann-Whitney rho_run, {res[0][0]} vs {D}: p={mannwhitney_p(a, b):.2g}")
