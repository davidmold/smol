#!/usr/bin/env python3
"""pw_compete2.py rebuild | setup | score   (output in ~/smol/runs/pwc2/)

Second round of the password contests (findings 37), after a review found two problems and
proposed a third contest:
- In 5 of 50 transitions a parasite-pool sequence was identical to a competing gene, adding
  69-155 extra copies of one side and biasing pB/pA in the direction predicted for escape.
  Here every pool sequence identical to any competitor in that transition is removed and the
  rest rescaled to 1,000 (removals are listed).
- 5 accepted transitions end at 2M, where persistence of the new password could not be checked;
  they are flagged ("unverified") and results are given with and without them.
- Reciprocal swap: the OLD password put on the NEW background (Mr, Gr), against the unmodified
  new pair. With the forward hybrid (new password on the old background) this separates the
  password's own effect from that of the background (hitchhiking predicts: the new pair keeps
  its advantage whichever password it carries; the forward hybrid is ~neutral).

Swaps: first the exact substitution used in round 1 (label in the machine; complement tag and the
label's codons in the gene, each found exactly once, "exact"). Evolved genes rarely carry the exact complement of the label (the label test aligns, so a partial
match is enough); then only the label's codons are swapped, plus the tag where it is present
exactly ("codons"). Kept only if the swapped pair closes a loop.

rebuild  reads pwc/transitions.json, writes pwc2/transitions.json
setup    contests, each of a FOCAL pair against an OTHER pair, 500 + 500 each, mutation off, 30k
         steps, 6 assay seeds, conditions p0 / pB / pA:
           main    focal (Mn, Gn)  vs other (Mo, Go)
           fwd     focal (Mh, Gh)  vs other (Mo, Go)    new password, old background
           recip   focal (Mn, Gn)  vs other (Mr, Gr)    new background: new vs old password
score    focal share of genes at 30k per transition (mean of seeds); per contest and condition the
         mean with a bootstrap 95% CI, sign counts, Wilcoxon; paired pB - p0 and pA - pB.
"""
import sys, os, re, json, random, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
H = os.path.expanduser("~"); R = f"{H}/smol/runs"; IN = f"{R}/pwc"; OUT = f"{R}/pwc2"
B = f"{H}/smol-src/target/release"
N, NP, STEPS, SEEDS = 500, 1000, 30000, (1, 2, 3, 4, 5, 6)

def lib():
    src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "pw_compete.py")).read()
    src = src.rsplit('\n{"find"', 1)[0]; g = {"__file__": __file__, "__name__": "lib"}; exec(src, g); return g

def swap(M, G, lfrom, lto, L):
    comp, enc, closes = L["comp"], L["enc"], L["closes"]
    Ms = M.replace(f"^F${lfrom}", f"^F${lto}", 1)
    if G.count(comp(lfrom)) == 1 and G.count(enc(lfrom)) == 1:
        g = G.replace(comp(lfrom), comp(lto)).replace(enc(lfrom), enc(lto))
        if closes(Ms, g): return Ms, g, "exact"
    # Evolved genes rarely carry the exact complement of the label: the label test aligns, so a
    # partial match finds the start. Then swap the label's codons only (so translation gives Ms)
    # and swap the tag too where it is present exactly.
    if G.count(enc(lfrom)) == 1:
        g = G.replace(enc(lfrom), enc(lto))
        if G.count(comp(lfrom)) == 1: g = g.replace(comp(lfrom), comp(lto))
        if closes(Ms, g): return Ms, g, "codons"
    return None, None, None

def rebuild():
    L = lib(); tr = json.load(open(f"{IN}/transitions.json")); os.makedirs(OUT, exist_ok=True)
    stats = collections.Counter()
    for x in tr:
        x["Mh"], x["Gh"], x["fwd_method"] = swap(x["Mo"], x["Go"], x["old_pw"], x["new_pw"], L)
        x["Mr"], x["Gr"], x["recip_method"] = swap(x["Mn"], x["Gn"], x["new_pw"], x["old_pw"], L)
        x["unverified"] = x["t2"] >= 2000000
        comps = {x[k] for k in ("Mo", "Go", "Mn", "Gn", "Mh", "Gh", "Mr", "Gr") if x.get(k)}
        for cond in ("pB", "pA"):
            rm = [(c, s) for c, s in x[cond] if s in comps]
            x[cond] = [(c, s) for c, s in x[cond] if s not in comps]
            x[cond + "_removed"] = len(rm)
            if rm: stats[f"{cond} pools with competitor sequences removed"] += 1
        stats[f"fwd {x['fwd_method']}"] += 1; stats[f"recip {x['recip_method']}"] += 1
    json.dump(tr, open(f"{OUT}/transitions.json", "w"), indent=1)
    print(len(tr), "transitions;", dict(stats), "; unverified persistence:", sum(x["unverified"] for x in tr))

def contests(x):
    c = {"main": ((x["Mn"], x["Gn"]), (x["Mo"], x["Go"]))}
    if x.get("Gh"): c["fwd"] = ((x["Mh"], x["Gh"]), (x["Mo"], x["Go"]))
    if x.get("Gr"): c["recip"] = ((x["Mn"], x["Gn"]), (x["Mr"], x["Gr"]))
    return c

def setup():
    tr = json.load(open(f"{OUT}/transitions.json")); jobs = open(f"{OUT}/jobs.txt", "w"); n = 0
    for k, x in enumerate(tr):
        base = [l for l in open(f"{R}/{x['set']}/{x['run']}/run.conf").read().splitlines()
                if l.split() and l.split()[0] not in ("RANDSEED", "NSTEPS", "MUTATE", "AGENT", "GRIDPOS", "CKPTEVERY",
                                                      "KEEPCKPT", "CKPTKEEPEVERY", "POPDYEVERY", "COPYEVERY", "SNAPEVERY")]
        for cname, ((Mf, Gf), (Mx, Gx)) in contests(x).items():
            for cond in ("p0", "pB", "pA"):
                if cond != "p0" and not x[cond]: continue
                for seed in SEEDS:
                    d = f"{OUT}/t{k:03d}/{cname}/{cond}/s{seed}"; os.makedirs(d, exist_ok=True)
                    mix = [(Mx, N), (Gx, N), (Mf, N), (Gf, N)]
                    if cond != "p0":
                        tot = sum(c for c, _ in x[cond]); mix += [(s, round(NP * c / tot)) for c, s in x[cond]]
                    rng = random.Random(seed * 1000 + 7)
                    cells = rng.sample([(i, j) for i in range(125) for j in range(100)], sum(c for _, c in mix))
                    lines = [f"RANDSEED {seed}", f"NSTEPS {STEPS}", "MUTATE 0"] + base + \
                            ["POPDYEVERY 1000", "SNAPEVERY 100000000", "CKPTEVERY 1000000"]
                    i = 0
                    for s, c in mix:
                        if c <= 0: continue
                        lines.append(f"AGENT {s} {c} Q")
                        for _ in range(c):
                            a, b = cells[i]; i += 1; lines.append(f"GRIDPOS {a} {b}")
                    open(f"{d}/run.conf", "w").write("\n".join(lines) + "\n")
                    jobs.write(f"cd {d} && {B}/smol run run.conf > log.txt 2>&1; rm -f copies.tsv grid.snap\n"); n += 1
    jobs.close(); print(n, "jobs")

def score():
    try:
        from scipy import stats
    except ImportError:
        py = f"{H}/smenv/bin/python"; os.execv(py, [py] + sys.argv)
    import numpy as np, statistics as S
    L = lib(); count = L["count"]; rng = np.random.default_rng(1)
    tr = json.load(open(f"{OUT}/transitions.json")); res = collections.defaultdict(dict); missing = 0; ext = 0
    for k, x in enumerate(tr):
        for cname, ((Mf, Gf), (Mx, Gx)) in contests(x).items():
            for cond in ("p0", "pB", "pA"):
                vals = []
                for seed in SEEDS:
                    d = f"{OUT}/t{k:03d}/{cname}/{cond}/s{seed}"
                    if not os.path.exists(f"{d}/run.conf"): continue
                    c = count(d) if os.path.exists(f"{d}/popdy001.dat") else None
                    if c is None:
                        log = open(f"{d}/log.txt", errors="replace").read() if os.path.exists(f"{d}/log.txt") else ""
                        if "EXTINCT" in log: ext += 1
                        else: missing += 1
                        continue
                    g = c[Gf] + c[Gx]
                    if g: vals.append(c[Gf] / g)
                if vals: res[(k, cname)][cond] = S.mean(vals)
    print(f"missing assay runs (technical): {missing}; extinct: {ext}")
    def summ(v, lab, ref=0.0):
        v = np.array(v)
        if len(v) < 3: return
        ci = stats.bootstrap((v,), np.mean, n_resamples=9999, method="percentile", random_state=rng).confidence_interval
        print(f"  {lab:34} n={len(v):2}  mean {v.mean():+.3f} [95% CI {ci.low:+.3f}, {ci.high:+.3f}]  "
              f"> {ref:g} in {int((v > ref).sum())}/{len(v)}  Wilcoxon p={stats.wilcoxon(v - ref).pvalue:.2g}")
    for subset in ("all", "verified"):
        for cname, lab in (("main", "new pair vs old pair"), ("fwd", "new password on OLD background vs old pair"),
                           ("recip", "new pair vs OLD password on NEW background")):
            rows = {k: v for (k, c), v in res.items() if c == cname and "p0" in v and (subset == "all" or not tr[k]["unverified"])}
            if len(rows) < 3: continue
            print(f"\n== {cname} ({lab}), {subset} transitions: focal share of genes at {STEPS // 1000}k")
            for cd in ("p0", "pB", "pA"):
                summ([v[cd] - .5 for v in rows.values() if cd in v], f"{cd} share - 0.5")
            summ([v["pB"] - v["p0"] for v in rows.values() if "pB" in v], "pB - p0 (vs earlier parasites)")
            summ([v["pA"] - v["pB"] for v in rows.values() if "pA" in v and "pB" in v], "pA - pB (later parasites)")
    json.dump({f"{k}:{c}": v for (k, c), v in res.items()}, open(f"{OUT}/scores.json", "w"), indent=1)

{"rebuild": rebuild, "setup": setup, "score": score}[sys.argv[1]]()
