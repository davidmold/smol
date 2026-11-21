#!/usr/bin/env python3
"""classify_gm.py SRC_DIR STEP: classify.py for gene-and-machine worlds (STRANDCOPY + TRANSLATE).

As classify.py, species are classed by what they did in the window after STEP (copies.tsv):
R made an exact copy (of a partner, or itself by translation); species copied exactly but never
copying are then split:
  G  a working gene: one of the world's three commonest machines, bound to it (smol-react --strand
     --translate 1 --stop, mutation off), copies it and translates it into a machine that, bound
     to it again, copies it AND makes itself (a closed gene-machine loop);
  P  anything else that machines copy: a gene-parasite (copied, but encodes no working machine).
N  neither. Writes classes.tsv (species, class, molecules at STEP), so swap.py ... P swaps only
gene-parasites, and prints the molecule shares."""
import sys, gzip, collections, subprocess, tempfile, os, re
d, step = sys.argv[1], int(sys.argv[2])
B = os.path.expanduser("~/smol-src/target/release")
mols = collections.Counter()
with gzip.open(f"{d}/checkpoint_{step}.txt.gz", "rt") as f:
    for l in f:
        if l.startswith("MOLS"): break
    for l in f:
        p = l.split()
        if len(p) > 3: mols[int(p[3])] += 1
seq = {}
for l in open(f"{d}/species.tsv"):
    p = l.rstrip("\n").split("\t")
    if p[0].isdigit() and len(p) >= 5: seq.setdefault(int(p[0]), p[-1])
copier, copied = set(), set()
for l in open(f"{d}/copies.tsv"):
    p = l.split("\t")
    if len(p) < 5 or not p[0].isdigit(): continue
    a, pa, c = int(p[1]), int(p[2]), int(p[3])
    if c == pa or c == a:
        copier.add(a); copied.add(c)

def react(m, g):
    """product sequences when machine m is bound to string g (m initiates)."""
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write(f"R M {m}\nR G {g}\n"); n = f.name
    out = subprocess.run([f"{B}/smol-react", n, "5000", "--strand", "--translate", "1", "--stop"],
                         capture_output=True, text=True).stdout
    os.unlink(n)
    block = re.search(r"^M\s+initiates with G.*?(?=^\S|\Z)", out, re.M | re.S)
    if not block: return []
    named = {"G": g, "M": m}
    res = []
    for l in block.group(0).splitlines():
        mm = re.match(r"\s+-> (?:product|passive made) (\S+)\s+(\S+)", l)
        if mm: res.append(named.get(mm.group(1), mm.group(2)))
    return res

_cache = {}
def working_gene(g, top):
    """g is a working gene if one of the world's three commonest machines, bound to it, copies it and
    translates it into a machine m2 that, bound to g, again copies g and makes m2."""
    if (g, top) in _cache: return _cache[(g, top)]
    ok = False
    if g and top:
        p1 = react(top, g)
        if g in p1:
            for m2 in [x for x in p1 if x != g]:
                p2 = react(m2, g)
                if g in p2 and m2 in p2: ok = True; break
    _cache[(g, top)] = ok
    return ok

reps = sorted([s for s in mols if s in copier and seq.get(s)], key=lambda s: -mols[s])
tops = [seq[s] for s in reps[:3]]   # the three commonest machines
cls = {}
for s in mols:
    if s in copier: cls[s] = "R"
    elif s in copied: cls[s] = "G" if any(working_gene(seq.get(s, ""), t) for t in tops) else "P"
    else: cls[s] = "N"
with open(f"{d}/classes.tsv", "w") as o:
    for s, n in mols.most_common():
        o.write(f"{s}\t{cls[s]}\t{n}\n")
tot = sum(mols.values()); by = collections.Counter()
for s, n in mols.items(): by[cls[s]] += n
print(d, tot, " ".join(f"{k}={by[k]} ({by[k] / tot:.0%})" for k in "RGPN"))
