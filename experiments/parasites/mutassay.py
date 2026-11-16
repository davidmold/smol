#!/usr/bin/env python3
"""mutassay.py paper|smol PATH... > table

Mutant-discrimination assay. For the top replicator of a world (the most abundant species that
copies itself in the reaction probe), make every single-symbol substitution and deletion, and
react the replicator with each mutant (smol-react --hub, mutation off). A mutant is
*functional* if it still copies itself, otherwise it is a would-be parasite. Reported:
  nf      share of mutants that are would-be parasites
  bindP   mean P(bind) of the replicator with would-be parasites
  bindF   the same with functional mutants
  copyP   of would-be parasites, binding-weighted share the replicator copies exactly
  copyF   the same for functional mutants (kin acceptance)
  disc    1 - copyP/copyF  (0: treats near-self parasites like kin; 1: never copies them)
  expP    per encounter with a random would-be parasite: P(bind) x copied (the parasite's success)

  paper RUN_DIR...   published runs (out1_STEP.conf every 100k)
  smol  RUN_DIR...   smol runs (checkpoint_STEP.txt.gz every 100k)
"""
import sys, os, re, glob, gzip, collections, subprocess, tempfile, statistics

REACT = os.path.expanduser("~/smol-src/target/release/smol-react")
KEY = "ABC$DEF%GH^IJK?LMN}OPQ>RST=UVWXYZ"
mode, paths = sys.argv[1], sys.argv[2:]

def react(seqs, flag, steps=5000):
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        for i, s in enumerate(seqs):
            f.write(f"R m{i} {s}\n")
        name = f.name
    out = subprocess.run([REACT, name, str(steps), flag], capture_output=True, text=True).stdout
    os.unlink(name)
    res = {}; cur = None
    for l in out.splitlines():
        m = re.match(r"(m\d+)\s+initiates with (m\d+)\s+P\(bind\)=([\d.]+)", l)
        if m:
            cur = (int(m.group(1)[1:]), int(m.group(2)[1:])); res[cur] = [float(m.group(3)), set()]
            continue
        m = re.match(r"\s+-> (?:product|passive made) (m\d+)", l)
        if m and cur:
            res[cur][1].add(int(m.group(1)[1:]))
    return res

def top_replicator(counts):
    top = [s for s, _ in counts.most_common(12)]
    res = react(top, "--self")
    for i, s in enumerate(top):
        if i in res.get((i, i), [0, set()])[1]:
            return s, counts[s] / sum(counts.values())
    return None, 0

def paper_snapshots(d):
    for f in sorted(glob.glob(f"{d}/out1_*.conf"), key=lambda p: int(re.search(r"_(\d+)\.", p).group(1))):
        step = int(re.search(r"_(\d+)\.", f).group(1))
        if step == 0 or step % 100000: continue
        c = collections.Counter()
        for l in open(f, errors="replace"):
            if l.startswith("AGENT"):
                c[l.split()[1]] += 1
            elif l.startswith("###"):
                p = l.split()
                if len(p) >= 3: c[p[2]] += 1
        yield step, c

def smol_snapshots(d):
    for f in sorted(glob.glob(f"{d}/checkpoint_*00000.txt.gz"), key=lambda p: int(re.search(r"_(\d+)\.", p).group(1))):
        step = int(re.search(r"_(\d+)\.", f).group(1))
        if step == 0: continue
        c = collections.Counter()
        with gzip.open(f, "rt") as fh:
            for l in fh:
                if l.startswith("MOLS"): break
            for l in fh:
                p = l.split()
                if p[5] == "F": c[p[4].replace(".", "")] += 1
        yield step, c

def assay(S):
    muts = set()
    for i in range(len(S)):
        for k in KEY:
            if k != S[i]: muts.add(S[:i] + k + S[i + 1:])
        muts.add(S[:i] + S[i + 1:])
    muts.discard(S); muts.discard("")
    muts = sorted(muts)
    res = react([S] + muts, "--hub")
    nP = nF = 0; bP = bF = 0.0; cP = cF = 0.0
    for j in range(1, len(muts) + 1):
        func = j in res.get((j, j), [0, set()])[1]
        p, prods = res.get((0, j), [0.0, set()])
        cop = j in prods
        if func: nF += 1; bF += p; cF += p * cop
        else:    nP += 1; bP += p; cP += p * cop
    copyP = cP / bP if bP else float("nan"); copyF = cF / bF if bF else float("nan")
    disc = 1 - copyP / copyF if copyF and copyF == copyF else float("nan")
    return dict(n=len(muts), nf=nP / len(muts), expP=cP / max(nP, 1), bindP=bP / max(nP, 1), bindF=bF / max(nF, 1),
                copyP=copyP, copyF=copyF, disc=disc)

rows = []
for d in paths:
    snaps = paper_snapshots(d) if mode == "paper" else smol_snapshots(d)
    label = "/".join(d.rstrip("/").split("/")[-2:])
    for step, c in snaps:
        S, share = top_replicator(c)
        if not S: 
            print(f"{label:16} {step:8d} no self-replicator in the top 12", flush=True); continue
        a = assay(S)
        rows.append((label, step, a))
        print(f"{label:16} {step:8d} len {len(S):3d} share {share:.2f} q {S.count('?'):2d}  nf {a['nf']:.2f} bindP {a['bindP']:.2f} bindF {a['bindF']:.2f}"
              f" copyP {a['copyP']:.2f} copyF {a['copyF']:.2f} disc {a['disc']:.2f} expP {a['expP']:.2f}  {S}", flush=True)

print("\n== medians by step")
by = collections.defaultdict(list)
for label, step, a in rows:
    by[(label.split("/")[0] if mode == "smol" else "paper", step)].append(a)
for k in sorted(by):
    v = by[k]; med = lambda key: statistics.median([x[key] for x in v if x[key] == x[key]] or [float("nan")])
    print(f"{k[0]:6} {k[1]:8d} n={len(v):2d}  nf {med('nf'):.2f} bindP {med('bindP'):.2f} bindF {med('bindF'):.2f} copyP {med('copyP'):.2f} copyF {med('copyF'):.2f} disc {med('disc'):.2f} expP {med('expP'):.2f}")
