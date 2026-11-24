#!/usr/bin/env python3
"""mech2.py SET [SET...] [-j N] -> direct reciprocal-adaptation assay on a set's classification
windows (~/smol/runs/ts_SET/sN/wTk: species.tsv + classes.tsv from classify_gm.py).

Replaces gmt_mech.py's thin probe (top machine, 5 parasites, two times). For every run and every
pair of window times (A, B), all 10 x 10:
  machines of A: the commonest machine species (class R) up to 80% of machine molecules (max 8)
  parasites of B: the commonest gene-parasites (P) up to 80% of parasite molecules (max 15)
  genes of B: the commonest working genes (G) up to 80% of gene molecules (max 5)
Each machine is reacted with each target in a fixed standard setting (smol-react --strand
--translate 1 --stop --hub, plus GMFLAGS, mutation off, 2 trials), abundance-weighted on both
sides:
  infect[A,B]  share of binding-weighted parasite encounters in which the parasite was copied
               (copying conditional on binding)
  gcopy[A,B]   the same for working genes: can A's machines still copy B's genes?
  bind[A,B]    mean P(bind) per encounter, parasites and genes separately: a defence can also
               work by binding less, which copying-given-binding would miss
Coverage: the share of each class's molecules the chosen species actually cover is recorded
(below 80% when the species cap binds) and reported.
Per run:
  host trend      for each parasite time B, Spearman of infect[A,B] with A; averaged over B.
                  Negative: later machines resist a fixed set of parasites better.
  parasite trend  for each machine time A, Spearman of infect[A,B] with B; averaged over A.
                  Positive: later parasites exploit a fixed set of machines better.
  gene trend      as host trend on gcopy. If later machines also lose the ability to copy
                  earlier GENES, resisting earlier parasites is a change of compatibility (a new
                  password), not a defence aimed at parasites; specificity = gcopy - infect.
An arms race in the strict sense needs a negative host trend AND a positive parasite trend, with
the host trend not explained by lost compatibility (specificity rising with A).
Writes ts_SET/mech2.json (all matrices) and prints a summary.
"""
import sys, os
try:
    import scipy  # noqa
except ImportError:
    py = os.path.expanduser("~/smenv/bin/python")
    if os.path.exists(py) and sys.executable != py: os.execv(py, [py] + sys.argv)
    raise
import re, glob, json, subprocess, tempfile, statistics as S
from multiprocessing import Pool
import numpy as np
from scipy import stats
R = os.path.expanduser("~/smol/runs"); B = os.path.expanduser("~/smol-src/target/release")
FLAGS = os.environ.get("GMFLAGS", "--qsearch").split()
CAPS = [int(x) for x in os.environ.get("MECH_CAPS", "8,15,5").split(",")]  # species caps: machines, parasites, genes

def top(lst, share, cap):
    """commonest species up to `share` of the class's molecules or `cap` species; returns them and
    the share of the class's molecules they actually cover (below `share` when the cap binds)."""
    lst = sorted(lst, reverse=True); tot = sum(c for c, _ in lst); out = []; acc = 0
    for c, s in lst:
        if acc >= share * tot or len(out) >= cap: break
        out.append((c, s)); acc += c
    return out, (acc / tot if tot else float("nan"))

def window(w):
    seq = {}
    for l in open(f"{w}/species.tsv"):
        p = l.rstrip("\n").split("\t")
        if p[0].isdigit() and len(p) >= 5: seq.setdefault(int(p[0]), p[-1])
    by = {"R": [], "P": [], "G": []}
    for l in open(f"{w}/classes.tsv"):
        a, b, c = l.split()
        if b in by and int(a) in seq: by[b].append((int(c), seq[int(a)]))
    (m, cm), (p, cp), (g, cg) = top(by["R"], .8, CAPS[0]), top(by["P"], .8, CAPS[1]), top(by["G"], .8, CAPS[2])
    return m, p, g, (cm, cp, cg)

def react(m, targets, trial):
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write(f"R M {m}\n" + "".join(f"R t{i} {s}\n" for i, s in enumerate(targets))); fn = f.name
    out = subprocess.run([f"{B}/smol-react", fn, "5000", "--strand", "--translate", "1", "--stop", "--hub",
                          "--seed", str(trial)] + FLAGS, capture_output=True, text=True).stdout
    os.unlink(fn); res = []
    for i in range(len(targets)):
        mm = re.search(rf"^M\s+initiates with t{i}\s+P\(bind\)=([\d.]+).*?(?=^\S|\Z)", out, re.M | re.S)
        if not mm: res.append((0.0, False)); continue
        res.append((float(mm.group(1)), bool(re.search(rf"-> (?:product|passive made) t{i}\b", mm.group(0)))))
    return res

def one_run(args):
    ts, r = args
    ws = sorted(glob.glob(f"{ts}/{r}/w*k"), key=lambda x: int(re.search(r"w(\d+)k", x).group(1)))
    W = {}
    for w in ws:
        if os.path.exists(f"{w}/classes.tsv"): W[int(re.search(r"w(\d+)k", w).group(1)) * 1000] = window(w)
    T = sorted(W); inf = {}; gc = {}; bd = {}; bg = {}
    for a in T:
        for b in T:
            ms, ps, gs = W[a][0], W[b][1], W[b][2]
            if not ms: continue
            targets = [s for _, s in ps] + [s for _, s in gs]; wts = [c for c, _ in ps] + [c for c, _ in gs]
            num = [0.0, 0.0]; den = [0.0, 0.0]; enc = [0.0, 0.0]
            for cm, m in ms:
                for trial in (1, 2):
                    for i, (pb, cp) in enumerate(react(m, targets, trial)):
                        k = 0 if i < len(ps) else 1; w = cm * wts[i] * pb
                        enc[k] += cm * wts[i]; den[k] += w; num[k] += w * cp
            if den[0]: inf[(a, b)] = num[0] / den[0]
            if den[1]: gc[(a, b)] = num[1] / den[1]
            if enc[0]: bd[(a, b)] = den[0] / enc[0]
            if enc[1]: bg[(a, b)] = den[1] / enc[1]
    cov = {t: W[t][3] for t in T}
    return r, T, inf, gc, bd, bg, cov

def trend(M, T, along):
    """mean over the fixed axis of Spearman(M, varying time)."""
    rs = []
    for fixed in T:
        pts = [(x, M[(x, fixed)] if along == "A" else M[(fixed, x)]) for x in T
               if ((x, fixed) if along == "A" else (fixed, x)) in M]
        if len(pts) >= 4 and len({v for _, v in pts}) > 1:
            rs.append(stats.spearmanr([p[0] for p in pts], [p[1] for p in pts]).statistic)
    return float(np.mean(rs)) if rs else float("nan")

def summarise(name, res):
    rows = {}
    for r, T, inf, gc, bd, bg, cov in res:
        spec = {k: gc[k] - inf[k] for k in inf if k in gc}
        rows[r] = {"host": trend(inf, T, "A"), "para": trend(inf, T, "B"), "gene": trend(gc, T, "A"),
                   "spec": trend(spec, T, "A"), "hbind": trend(bd, T, "A"), "gbind": trend(bg, T, "A"),
                   "covM": float(np.median([c[0] for c in cov.values()])), "covP": float(np.median([c[1] for c in cov.values()])),
                   "covG": float(np.median([c[2] for c in cov.values()])),
                   "contemp": float(np.mean([inf[(t, t)] for t in T if (t, t) in inf])) if T else float("nan")}
    print(f"\n== {name}: {len(rows)} runs")
    for key, lab, esc in (("host", "host trend (later machines vs fixed parasites)", -1),
                          ("para", "parasite trend (later parasites vs fixed machines)", +1),
                          ("gene", "gene trend (later machines vs fixed genes)", -1),
                          ("spec", "specificity trend (gene minus parasite copying)", +1),
                          ("hbind", "binding trend (later machines binding fixed parasites)", -1),
                          ("gbind", "binding trend (later machines binding fixed genes)", -1)):
        v = [x[key] for x in rows.values() if not np.isnan(x[key])]
        if len(v) < 3: continue
        k = sum(1 for x in v if x * esc > 0); nz = [x for x in v if x != 0]
        print(f"  {lab:52} median {S.median(v):+.3f}  arms-race sign {k}/{len(v)}  "
              f"sign p={stats.binomtest(k, len(v)).pvalue:.2g}  Wilcoxon p={stats.wilcoxon(nz).pvalue:.2g}")
    both = sum(1 for x in rows.values() if x["host"] < 0 and x["para"] > 0)
    print(f"  runs with host trend < 0 and parasite trend > 0: {both}/{len(rows)}")
    print(f"  contemporary exploitation infect[t,t], median over runs: {S.median([x['contemp'] for x in rows.values()]):.2f}")
    print("  abundance actually covered (median over runs of each run's median window): machines "
          + "  ".join(f"{S.median([x[k] for x in rows.values()]):.2f}" for k in ("covM",))
          + f"  parasites {S.median([x['covP'] for x in rows.values()]):.2f}  genes {S.median([x['covG'] for x in rows.values()]):.2f}")
    for r, x in rows.items():
        print(f"    {r}: " + "  ".join(f"{k} {x[k]:+.2f}" for k in ("host", "para", "gene", "spec", "hbind", "contemp", "covM", "covP")))
    return rows

if __name__ == "__main__":
    args = sys.argv[1:]; j = 4; allrows = []
    if "-j" in args: i = args.index("-j"); j = int(args[i + 1]); del args[i:i + 2]
    for a in args:
        name = a[3:] if a.startswith("ts_") else a; ts = f"{R}/ts_{name}"
        runs = sorted((os.path.basename(x) for x in glob.glob(f"{ts}/s*")), key=lambda x: int(x[1:]))
        with Pool(j) as pool: res = pool.map(one_run, [(ts, r) for r in runs])
        J = lambda M: {f"{a_}:{b_}": v for (a_, b_), v in M.items()}
        json.dump({r: {"times": T, "infect": J(inf), "gcopy": J(gc), "bind_p": J(bd), "bind_g": J(bg),
                       "coverage": {str(t): c for t, c in cov.items()}} for r, T, inf, gc, bd, bg, cov in res},
                  open(f"{ts}/mech2" + os.environ.get("MECH_TAG", "") + ".json", "w"))
        allrows.append((name, summarise(name, res)))
    for key in ("host", "para", "spec"):
        if len(allrows) < 2: break
        a = [x[key] for x in allrows[0][1].values() if not np.isnan(x[key])]
        for n2, rows in allrows[1:]:
            b = [x[key] for x in rows.values() if not np.isnan(x[key])]
            if len(a) >= 3 and len(b) >= 3:
                u = stats.mannwhitneyu(a, b, alternative="two-sided")
                print(f"{key} trend, {allrows[0][0]} vs {n2}: medians {S.median(a):+.3f} vs {S.median(b):+.3f}  "
                      f"Mann-Whitney p={u.pvalue:.2g}  rank-biserial {1 - 2 * u.statistic / (len(a) * len(b)):+.2f}")
