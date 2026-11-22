#!/usr/bin/env python3
"""gmt_mech.py TS_DIR -> mechanism of the tail-check arms race.
For each run and for times 200k, 1M, 2M (windows classified by classify_gm.py): the share of
gene-parasite molecules that carry the tail tag the machine searches for (complement of the
machine's second label test, found by `?` + letters after the first test), and a 2x2 probe:
top machine of time A against the 5 commonest gene-parasites of time B (A, B in {200k, 2M}),
abundance- and binding-weighted share copied (smol-react --qsearch, 3 trials)."""
import sys, os, re, glob, subprocess, tempfile, statistics as S, collections
TS = sys.argv[1]; B = os.path.expanduser("~/smol-src/target/release")
comp = lambda s: "".join(chr((ord(c) - 65 + 13) % 26 + 65) for c in s)
def window(r, t):
    w = f"{TS}/{r}/w{t // 1000}k"; seq = {}
    for l in open(f"{w}/species.tsv"):
        p = l.rstrip("\n").split("\t")
        if p[0].isdigit() and len(p) >= 5: seq.setdefault(int(p[0]), p[-1])
    cl = [(int(a), b, int(c)) for a, b, c in (l.split() for l in open(f"{w}/classes.tsv"))]
    reps = sorted([(c, a) for a, b, c in cl if b == "R" and a in seq], reverse=True)
    ps = sorted([(c, a) for a, b, c in cl if b == "P" and a in seq], reverse=True)
    return seq[reps[0][1]] if reps else None, [(c, seq[a]) for c, a in ps]
def copied(M, ps):
    ps = ps[:5]; tc = tt = 0
    for trial in (1, 2, 3):
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
            f.write(f"R M {M}\n" + "".join(f"R p{i} {s}\n" for i, (n, s) in enumerate(ps))); fn = f.name
        out = subprocess.run([f"{B}/smol-react", fn, "5000", "--strand", "--translate", "1", "--stop", "--hub", "--seed", str(trial)] + os.environ.get("MECHFLAGS", "--qsearch").split(), capture_output=True, text=True).stdout
        os.unlink(fn)
        for i, (n, s) in enumerate(ps):
            m = re.search(rf"^M\s+initiates with p{i}\s+P\(bind\)=([\d.]+).*?(?=^\S|\Z)", out, re.M | re.S)
            if not m: continue
            pb = float(m.group(1)); tt += n * pb
            if re.search(rf"-> (?:product|passive made) p{i}\b", m.group(0)): tc += n * pb
    return tc / tt if tt else float("nan")
tag = collections.defaultdict(list); probe = collections.defaultdict(list)
for r in sorted(os.path.basename(x) for x in glob.glob(f"{TS}/s*")):
    W = {t: window(r, t) for t in (200000, 1000000, 2000000)}
    for t, (M, ps) in W.items():
        labs = re.findall(r"\?([A-Z]{2,})", M or "")
        if len(labs) >= 2 and ps:
            tg = comp(labs[1]); n = sum(c for c, _ in ps)
            tag[t].append(sum(c for c, s in ps if tg in s) / n)
    for a in (200000, 2000000):
        for b in (200000, 2000000):
            if W[a][0] and W[b][1]: probe[(a, b)].append(copied(W[a][0], W[b][1]))
print("share of gene-parasite molecules carrying the machine's tail tag (median over runs):")
for t in (200000, 1000000, 2000000):
    if tag[t]: print(f"  {t // 1000:>5}k  {S.median(tag[t]):.2f}  (runs {len(tag[t])})")
print("top machine of time A vs gene-parasites of time B: share copied (median over runs)")
for (a, b), v in sorted(probe.items()):
    print(f"  machine {a // 1000:>5}k  parasites {b // 1000:>5}k   {S.median(v):.2f}   range {min(v):.2f}-{max(v):.2f}")
