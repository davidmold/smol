#!/usr/bin/env python3
"""pw_replay.py setup [--validate] | run CASE TYPE REP | score     (output in ~/smol/runs/pwr/)

Lineage replays (proposed by the second review of findings 37): does a rare founder leave a
lineage that establishes in the earlier, evolving world, and does that depend on its background
or on its translation-start label?

Cases: the transitions with all four reconstructed pairs (pwc2/transitions.json):
  ev_new  evolved background, new label   (Mn, Gn)
  ev_old  evolved background, old label   (Mr, Gr)   reciprocal swap
  an_new  ancestral background, new label (Mh, Gh)   forward swap
  an_old  ancestral background, old label (Mo, Go)
World: the run's own checkpoint at the earlier window t (the old pair dominant), MUTATION ON (the
run's own rate), resumed for 50k steps. Founders: 5 machine-gene pairs put at 5 matched sites
(machine at (x, y), gene at (x+1, y); the same sites for all four types in a replicate).

Tracking: each founder gene carries a synonymous marker, one codon in the middle of its coding
part changed to another codon for the same symbol (translation checked to give the same machine;
the pair checked to close a loop), so that the founder gene is a new species even where the same
sequence lives in the world. The same kind of marker is used for all four types. species.tsv
records every new species' template, so the founder's lineage is every species whose template
chain reaches the founder gene. Machines are not attributed (a translated machine may be a
species that already lives in the world); the lineage is scored by its gene-like strands (letters
only, at least 90% of the founder's length). Recorded at +10k, +25k, +50k: lineage gene-like
molecules, all lineage letter strands, and the share of lineage gene-like molecules that still
carry the founder label's codons (label kept) - so establishing with the original label can be
told from establishing after the label changes.
Established: >= 20 lineage gene-like molecules at +50k.

--validate: 1 case x 4 types x 2 reps with MUTATE 0 and 10k steps, to check the tracking (with
mutation off, the lineage must be exactly the marked founder gene).
"""
import sys, os, re, json, random, subprocess, tempfile, glob, collections, gzip, shutil
H = os.path.expanduser("~"); R = f"{H}/smol/runs"; B = f"{H}/smol-src/target/release"
T = f"{H}/smol-src/tools"; E = os.path.dirname(os.path.abspath(__file__)); OUT = f"{R}/pwr"
REPS, STEPS, NPAIR, EST = 20, 50000, 5, 20
TYPES = {"ev_new": ("Mn", "Gn"), "ev_old": ("Mr", "Gr"), "an_new": ("Mh", "Gh"), "an_old": ("Mo", "Go")}
LAB = {"ev_new": "new_pw", "ev_old": "old_pw", "an_new": "new_pw", "an_old": "old_pw"}
src = open(f"{E}/pw_compete.py").read().rsplit('\n{"find"', 1)[0]
L = {"__file__": f"{E}/pw_compete.py", "__name__": "lib"}; exec(src, L)
closes, enc, react = L["closes"], L["enc"], L["react"]

def table():
    t = {}
    for l in subprocess.run([f"{B}/smol-codon", "table"], capture_output=True, text=True).stdout.splitlines():
        p = l.split()
        if p: t[p[0]] = p[1:]
    return t
def decode(g): return subprocess.run([f"{B}/smol-codon", "decode", g], capture_output=True, text=True).stdout.strip()
def start(label, g):
    o = subprocess.run([f"{B}/smol-search"], input=f"{label} {g}\n", capture_output=True, text=True).stdout.split()
    return int(o[0])

def mark(M, G, label, tab):
    """G with one codon in the middle of its coding part swapped for a synonymous one."""
    p = start(label, G); base = decode(G[p:]); n = (len(G) - p) // 2
    for off in range(0, n // 2):
        for i in (n // 2 + off, n // 2 - off - 1):
            if i < 0 or i >= n: continue
            c = G[p + 2 * i:p + 2 * i + 2]
            sym = next((s for s, cs in tab.items() if c in cs), None)
            if sym is None: continue
            for alt in tab[sym]:
                if alt == c: continue
                g = G[:p + 2 * i] + alt + G[p + 2 * i + 2:]
                if decode(g[p:]) == base and start(label, g) == p and closes(M, g): return g, i
    return None, None

def setup():
    validate = "--validate" in sys.argv
    tr = json.load(open(f"{R}/pwc2/transitions.json")); tab = table(); os.makedirs(OUT, exist_ok=True)
    cases = {}
    for k, x in enumerate(tr):
        if not (x.get("Gh") and x.get("Gr")): continue
        known = set()
        for l in open(f"{R}/{x['set']}/{x['run']}/species.tsv"):
            p = l.rstrip("\n").split("\t")
            if len(p) >= 5: known.add(p[4])
        c = {"k": k, "set": x["set"], "run": x["run"], "t": x["t"], "old_pw": x["old_pw"], "new_pw": x["new_pw"], "types": {}}
        for ty, (mk, gk) in TYPES.items():
            M, G = x[mk], x[gk]; lab = x[LAB[ty]]
            g, i = mark(M, G, lab, tab)
            if g is None or g in known: c["types"][ty] = None; continue
            c["types"][ty] = {"M": M, "G": g, "G_unmarked": G, "label": lab, "marker_codon": i, "label_codons": enc(lab)}
        cases[k] = c
        print(f"t{k:03d}: " + "  ".join(f"{ty} {'ok' if v else 'NO MARKER'}" for ty, v in c["types"].items()))
        if validate: break
    json.dump(cases, open(f"{OUT}/cases{'_validate' if validate else ''}.json", "w"), indent=1)
    reps = 2 if validate else REPS
    with open(f"{OUT}/jobs{'_validate' if validate else ''}.txt", "w") as f:
        for k, c in cases.items():
            for ty, v in c["types"].items():
                if v is None: continue
                for r in range(1, reps + 1):
                    f.write(f"python3 {E}/pw_replay.py run {k} {ty} {r}{' --validate' if validate else ''}\n")
    print(sum(1 for _ in open(f"{OUT}/jobs{'_validate' if validate else ''}.txt")), "jobs")

def lineage_scan(d, founder, t0, times):
    sp = {}
    for l in open(f"{d}/species.tsv"):
        p = l.rstrip("\n").split("\t")
        if len(p) >= 5 and p[0].isdigit(): sp[int(p[0])] = (int(p[1]), int(p[3]), p[4])
    roots = {i for i, v in sp.items() if v[2] == founder}
    lin = set(roots); changed = True
    while changed:
        changed = False
        for i, (born, pp, s) in sp.items():
            if i not in lin and born >= t0 and pp in lin:
                lin.add(i); changed = True
    counts = {t: collections.Counter() for t in times}
    for l in open(f"{d}/popdy001.dat"):
        st, s, n = (int(v) for v in l.split(","))
        if st in counts and s in lin: counts[st][s] += n
    return sp, lin, counts

def run():
    k, ty, rep = sys.argv[2], sys.argv[3], int(sys.argv[4]); validate = "--validate" in sys.argv
    c = json.load(open(f"{OUT}/cases{'_validate' if validate else ''}.json"))[k]; v = c["types"][ty]
    d = f"{OUT}{'/validate' if validate else ''}/t{int(k):03d}/{ty}/r{rep}"
    if os.path.exists(f"{d}/result.json"): return
    shutil.rmtree(d, ignore_errors=True)
    t0 = c["t"]; steps = 10000 if validate else STEPS
    rng = random.Random(int(k) * 1000 + rep)
    sites = rng.sample([(x, y) for x in range(124) for y in range(100)], NPAIR)
    ops = []
    for x, y in sites: ops += ["--put", f"{x},{y}", v["M"], "--put", f"{x + 1},{y}", v["G"]]
    ops += ["--set", "NSTEPS", str(t0 + steps), "--set", "POPDYEVERY", "1000", "--set", "COPYEVERY", "100000000",
            "--set", "SPECIESFORGET", "1", "--set", "CKPTEVERY", "100000000"]
    if validate: ops += ["--set", "MUTATE", "0"]
    src = f"{R}/{c['set']}/{c['run']}"
    subprocess.run(["python3", f"{T}/ckpt_edit.py", src, f"checkpoint_{t0}.txt.gz", d, "--live-species", "--rng-seed", str(rep)] + ops,
                   check=True, capture_output=True)
    subprocess.run([f"{B}/smol", "resume", f"checkpoint_{t0}.txt.gz", "--reseed", str(rep)], cwd=d,
                   stdout=open(f"{d}/log.txt", "w"), stderr=subprocess.STDOUT)
    times = [t0 + s for s in ((2000, 5000, 10000) if validate else (10000, 25000, 50000))]
    sp, lin, counts = lineage_scan(d, v["G"], t0, times)
    L0 = len(v["G"]); res = {"case": int(k), "type": ty, "rep": rep, "sites": sites, "lineage_species": len(lin), "at": {}}
    for t in times:
        gl = {s: n for s, n in counts[t].items() if re.fullmatch(r"[A-Z]+", sp[s][2]) and len(sp[s][2]) >= 0.9 * L0}
        letters = sum(n for s, n in counts[t].items() if re.fullmatch(r"[A-Z]+", sp[s][2]))
        kept = sum(n for s, n in gl.items() if v["label_codons"] in sp[s][2])
        exact = counts[t].get(next(iter(i for i in lin if sp[i][2] == v["G"]), -1), 0)
        res["at"][t - t0] = {"genelike": sum(gl.values()), "letters": letters, "label_kept": kept, "exact_founder": exact,
                             "other_lineage": sum(n for s, n in counts[t].items() if sp[s][2] != v["G"])}
    last = res["at"][steps]; res["established"] = last["genelike"] >= EST
    json.dump(res, open(f"{d}/result.json", "w"))
    for f in os.listdir(d):
        if f not in ("result.json", "log.txt"): os.remove(f"{d}/{f}")

def score():
    try:
        from scipy import stats
    except ImportError:
        py = f"{H}/smenv/bin/python"; os.execv(py, [py] + sys.argv)
    import statistics as S
    validate = "--validate" in sys.argv; base = f"{OUT}/validate" if validate else OUT
    cases = json.load(open(f"{OUT}/cases{'_validate' if validate else ''}.json"))
    res = [json.load(open(f)) for f in glob.glob(f"{base}/t*/*/r*/result.json")]
    if validate:
        for r in sorted(res, key=lambda r: (r["type"], r["rep"])):
            print(r["type"], r["rep"], {t: (a["exact_founder"], a["other_lineage"], a["genelike"]) for t, a in r["at"].items()},
                  "lineage species", r["lineage_species"])
        return
    by = collections.defaultdict(list)
    for r in res: by[(r["case"], r["type"])].append(r)
    print(f"{len(res)} replays. Established = >= {EST} lineage gene-like molecules at +{STEPS // 1000}k (of {REPS} replicates)\n")
    print("case                     " + "".join(f"{ty:>16}" for ty in TYPES))
    rows = {}
    for k, c in sorted(cases.items(), key=lambda kv: int(kv[0])):
        k = int(k); row = {}
        for ty in TYPES:
            rs = by.get((k, ty), [])
            if rs:
                est = sum(r["established"] for r in rs); g = S.median(r["at"][str(STEPS)]["genelike"] for r in rs)
                kept = [r["at"][str(STEPS)]["label_kept"] / r["at"][str(STEPS)]["genelike"] for r in rs if r["at"][str(STEPS)]["genelike"]]
                row[ty] = (est / len(rs), g, S.median(kept) if kept else float("nan"), len(rs))
        rows[k] = row
        print(f"t{k:03d} {c['old_pw']:>5}->{c['new_pw']:<5} {c['set']:6} " + "".join(
            f"{row[ty][0]:>7.2f} ({row[ty][1]:>4.0f})" if ty in row else f"{'-':>16}" for ty in TYPES))
    print("\n(share of replicates established, and median lineage gene-like molecules at the end)")
    print("label kept among lineage gene-like molecules at the end (median over established and unestablished):")
    for k, row in rows.items():
        print(f"  t{k:03d} " + "  ".join(f"{ty} {row[ty][2]:.2f}" for ty in TYPES if ty in row))
    full = [k for k, row in rows.items() if all(ty in row for ty in TYPES)]
    if len(full) >= 3:
        print(f"\nAcross the {len(full)} cases with all four types (establishment share, paired by case):")
        for a, b, lab in (("ev_new", "an_new", "evolved vs ancestral background, new label"),
                          ("ev_old", "an_old", "evolved vs ancestral background, old label"),
                          ("ev_new", "ev_old", "new vs old label, evolved background"),
                          ("an_new", "an_old", "new vs old label, ancestral background")):
            dd = [rows[k][a][0] - rows[k][b][0] for k in full]
            nz = [v for v in dd if v != 0]
            p = stats.wilcoxon(nz).pvalue if len(nz) >= 3 else float("nan")
            print(f"  {lab:46} mean difference {sum(dd) / len(dd):+.2f}  positive in {sum(1 for v in dd if v > 0)}/{len(dd)}"
                  f"  zero in {sum(1 for v in dd if v == 0)}  Wilcoxon p={p:.2g}")

{"setup": setup, "run": run, "score": score}[sys.argv[1]]()
