#!/usr/bin/env python3
"""pw_compete.py find | setup | score   (output in ~/smol/runs/pwc/)

Do parasites select for password changes? (findings 36: machines keep changing their
translation-start password, and later parasites match later machines; what drives it is open.)

find   For every run of cf_gmt, cf_gmq and cf_gmd, the first change of the commonest machine's
       password between consecutive classification windows t -> t' (pw_scan.py) where the new
       password is still the commonest machine's at the next window (if there is one). For each:
         old pair (Mo, Go)   commonest machine at t, and the commonest working gene at t with
                             which it closes a loop (Mo bound to Go copies Go and makes Mo)
         new pair (Mn, Gn)   the same at t'
         hybrid (Mh, Gh)     Mo and Go with only the password changed to the new one (the label
                             in the machine, its complement tag in the gene, and the label's
                             codons in the gene); kept only if it closes a loop
         parasites before    gene-parasites (class P) at t, commonest up to 80% (max 10)
         parasites after     the same at t'
       Pairs that do not close a loop are dropped (with the reason). One transition per run, so
       transitions are independent replicates. Writes transitions.json.
setup  Contests in fresh worlds with that run's rules (its run.conf), mutation off, 30k steps:
         main    500 Mo + 500 Go + 500 Mn + 500 Gn
         hybrid  500 Mo + 500 Go + 500 Mh + 500 Gh   (the password change alone)
       each under three conditions: p0 no parasites, pB 1,000 parasites from before (split by
       abundance), pA 1,000 parasites from after; 4 assay seeds; random cells (same per seed).
score  Per transition and contest: the new side's share of genes and of machines (new / (new + old),
       exact sequences, mutation off) at 30k, mean over seeds; parasite molecules at 30k / 1,000;
       total genes at 30k. Across transitions (the replicates): the new gene share under each
       condition, and paired differences pB - p0 (escape advantage against old parasites) and
       pA - pB (erosion by new parasites), with Wilcoxon signed-rank tests.
Readings (from the review that proposed this):
  new side gains specifically against earlier parasites (pB > p0) and later parasites erode it
  (pA < pB): parasite-driven escape and counter-adaptation; new side wins equally without
  parasites (p0 > 0.5, pB ~ p0): a general reproductive advantage; no reproducible advantage:
  drift or unmeasured ecological context remains plausible.
"""
import sys, os, re, json, glob, random, subprocess, tempfile, collections
H = os.path.expanduser("~"); B = f"{H}/smol-src/target/release"; R = f"{H}/smol/runs"; OUT = f"{R}/pwc"
SETS = ("cf_gmt", "cf_gmq", "cf_gmd"); N, NP, STEPS, SEEDS = 500, 1000, 30000, (1, 2, 3, 4)
comp = lambda s: "".join(chr((ord(c) - 65 + 13) % 26 + 65) for c in s)

def react(m, g):
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write(f"R M {m}\nR G {g}\n"); fn = f.name
    out = subprocess.run([f"{B}/smol-react", fn, "5000", "--strand", "--translate", "1", "--stop", "--qsearch"],
                         capture_output=True, text=True).stdout
    os.unlink(fn)
    blk = re.search(r"^M\s+initiates with G.*?(?=^\S|\Z)", out, re.M | re.S)
    if not blk: return []
    named = {"G": g, "M": m}
    return [named.get(a, b) for a, b in re.findall(r"-> (?:product|passive made) (\S+)\s+(\S+)", blk.group(0))]
def closes(m, g):
    p = react(m, g); return g in p and m in p
def enc(s): return subprocess.run([f"{B}/smol-codon", "encode", s], capture_output=True, text=True).stdout.strip()
def label(m):
    x = re.search(r"\^F\$([A-Z]+)", m or ""); return x.group(1) if x else None

def pair(by, M):
    """M with the commonest working gene it closes a loop with (genes of class G, then any)."""
    for _, g in by["G"][:8]:
        if closes(M, g): return g
    return None

def find():
    sys.argv = ["x"]; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import pw_scan
    os.makedirs(OUT, exist_ok=True); tr = []; dropped = collections.Counter()
    for name in SETS:
        for r, out in pw_scan.runs(name):
            for i in range(len(out) - 1):
                (t, Mo, lo, _, bo), (t2, Mn, ln, _, bn) = out[i], out[i + 1]
                if not (lo and ln and lo != ln): continue
                if i + 2 < len(out) and out[i + 2][2] != ln: dropped["new password not kept"] += 1; continue
                Go, Gn = pair(bo, Mo), pair(bn, Mn)
                if not Go or not Gn: dropped["pair does not close a loop"] += 1; continue
                Mh = Mo.replace(f"^F${lo}", f"^F${ln}", 1); Gh = None
                if Go.count(comp(lo)) == 1 and Go.count(enc(lo)) == 1:
                    g = Go.replace(comp(lo), comp(ln)).replace(enc(lo), enc(ln))
                    if closes(Mh, g): Gh = g
                if Gh is None: dropped["hybrid not built (kept main contest)"] += 1
                def pool(by):
                    ps = by["P"]; tot = sum(c for c, _ in ps); acc = 0; sel = []
                    for c, s in ps:
                        if acc >= .8 * tot or len(sel) >= 10: break
                        sel.append((c, s)); acc += c
                    return sel
                tr.append({"set": name, "run": r, "t": t, "t2": t2, "old_pw": lo, "new_pw": ln,
                           "Mo": Mo, "Go": Go, "Mn": Mn, "Gn": Gn, "Mh": Mh, "Gh": Gh,
                           "pB": pool(bo), "pA": pool(bn)})
                break
    json.dump(tr, open(f"{OUT}/transitions.json", "w"), indent=1)
    print(f"{len(tr)} transitions ({collections.Counter(x['set'] for x in tr)}); hybrids {sum(1 for x in tr if x['Gh'])}")
    print("dropped:", dict(dropped))

def setup():
    tr = json.load(open(f"{OUT}/transitions.json")); jobs = open(f"{OUT}/jobs.txt", "w"); n = 0
    for k, x in enumerate(tr):
        base = [l for l in open(f"{R}/{x['set']}/{x['run']}/run.conf").read().splitlines()
                if l.split() and l.split()[0] not in ("RANDSEED", "NSTEPS", "MUTATE", "AGENT", "GRIDPOS", "CKPTEVERY",
                                                      "KEEPCKPT", "CKPTKEEPEVERY", "POPDYEVERY", "COPYEVERY", "SNAPEVERY")]
        contests = {"main": (x["Mn"], x["Gn"])}
        if x["Gh"]: contests["hybrid"] = (x["Mh"], x["Gh"])
        for cname, (M2, G2) in contests.items():
            for cond in ("p0", "pB", "pA"):
                for seed in SEEDS:
                    d = f"{OUT}/t{k:03d}/{cname}/{cond}/s{seed}"; os.makedirs(d, exist_ok=True)
                    mix = [(x["Mo"], N), (x["Go"], N), (M2, N), (G2, N)]
                    if cond != "p0":
                        ps = x[cond]; tot = sum(c for c, _ in ps)
                        mix += [(s, round(NP * c / tot)) for c, s in ps]
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

def count(d, step=STEPS):
    seq = {}
    for l in open(f"{d}/species.tsv"):
        p = l.rstrip("\n").split("\t")
        if p[0].isdigit() and len(p) >= 5: seq[int(p[0])] = p[-1]
    by = collections.Counter(); ok = False
    for l in open(f"{d}/popdy001.dat"):
        st, sp, c = (int(v) for v in l.split(","))
        if st == step: by[seq.get(sp, "?")] += c; ok = True
    return by if ok else None

def score():
    import statistics as S
    try:
        from scipy import stats
    except ImportError:
        py = f"{H}/smenv/bin/python"; os.execv(py, [py] + sys.argv)
    tr = json.load(open(f"{OUT}/transitions.json")); res = collections.defaultdict(dict); missing = 0
    for k, x in enumerate(tr):
        for cname, (M2, G2) in (("main", (x["Mn"], x["Gn"])), ("hybrid", (x["Mh"], x["Gh"]))):
            if not G2: continue
            for cond in ("p0", "pB", "pA"):
                gs, ms, ps, tot = [], [], [], []
                for seed in SEEDS:
                    d = f"{OUT}/t{k:03d}/{cname}/{cond}/s{seed}"
                    c = count(d) if os.path.exists(f"{d}/popdy001.dat") else None
                    if c is None:
                        log = open(f"{d}/log.txt", errors="replace").read() if os.path.exists(f"{d}/log.txt") else ""
                        if "EXTINCT" in log: gs.append(float("nan")); continue
                        missing += 1; continue
                    g = c[G2] + c[x["Go"]]; m = c[M2] + c[x["Mo"]]
                    if g: gs.append(c[G2] / g)
                    if m: ms.append(c[M2] / m)
                    pset = {s for _, s in x["pB"]} | {s for _, s in x["pA"]}
                    if cond != "p0": ps.append(sum(c[s] for s in pset) / NP)
                    tot.append(g)
                gs = [v for v in gs if v == v]
                res[(k, cname)][cond] = {"gene": S.mean(gs) if gs else float("nan"), "mach": S.mean(ms) if ms else float("nan"),
                                         "par": S.mean(ps) if ps else float("nan"), "genes": S.mean(tot) if tot else float("nan")}
    print(f"missing assay runs (technical): {missing}")
    for cname in ("main", "hybrid"):
        rows = {k: v for (k, c), v in res.items() if c == cname and all(cd in v for cd in ("p0", "pB", "pA"))}
        for setname in ("all",) + SETS:
            ks = [k for k in rows if setname == "all" or tr[k]["set"] == setname]
            if len(ks) < 3: continue
            print(f"\n== {cname} contest, {setname}: {len(ks)} transitions (new side's share of genes at {STEPS // 1000}k)")
            for cd in ("p0", "pB", "pA"):
                v = [rows[k][cd]["gene"] for k in ks if rows[k][cd]["gene"] == rows[k][cd]["gene"]]
                if len(v) >= 3:
                    print(f"  {cd}: median {S.median(v):.2f}  new side ahead in {sum(1 for a in v if a > .5)}/{len(v)}  "
                          f"Wilcoxon vs 0.5 p={stats.wilcoxon([a - .5 for a in v]).pvalue:.2g}")
            for a, b, lab in (("pB", "p0", "pB - p0 (advantage against earlier parasites)"),
                              ("pA", "pB", "pA - pB (erosion by later parasites)")):
                dd = [rows[k][a]["gene"] - rows[k][b]["gene"] for k in ks
                      if rows[k][a]["gene"] == rows[k][a]["gene"] and rows[k][b]["gene"] == rows[k][b]["gene"]]
                if len(dd) >= 3:
                    print(f"  {lab}: median {S.median(dd):+.3f}  positive in {sum(1 for v in dd if v > 0)}/{len(dd)}  "
                          f"Wilcoxon p={stats.wilcoxon(dd).pvalue:.2g}")
            for cd in ("pB", "pA"):
                v = [rows[k][cd]["par"] for k in ks if rows[k][cd]["par"] == rows[k][cd]["par"]]
                if v: print(f"  parasites at end / start, {cd}: median {S.median(v):.2f}")
    json.dump({f"{k}:{c}": v for (k, c), v in res.items()}, open(f"{OUT}/scores.json", "w"), indent=1)

{"find": find, "setup": setup, "score": score}[sys.argv[1]]()
