#!/usr/bin/env python3
"""What replaced the original pair? Ancestry plus the placed-copy network at the end of a run.

    lineage.py RUN_DIR NAME=SEQ [NAME=SEQ ...] [--window 20000] [--json OUT]

Lineage: every species is assigned to the founder it descends from by following its template
parent (species.tsv column 4, the passive at the product's reaction: in the copy-strand
chemistry the template) back to a seeded sequence. Named sequences are the founders.

Reported for the last WINDOW steps, over placed exact copies (product == template):
  pop     share of molecules at the end by lineage
  copies  share of placed exact copies by (maker lineage -> template lineage), with
          'self' = maker and template the same species
  cross   placed exact copies between the A and B lineages in each direction, and how many
          distinct species pairs carry them (1-2 = one partnership; many = interchangeable)
  roles   among species with >= 1% of molecules: lineage, length, has '=', placed copies made
          as maker, placed exact copies of it made by itself / by its own lineage / by others
"""
import argparse, collections, gzip, json, os, sys

def species_lines(p):
    """Lines of species.tsv and species.tsv.gz together (a run resumed from an edited checkpoint
    with GZIPLOGS has both: the copied list and the species born after the resume)."""
    import gzip as _g, os as _o
    for q in (p, p + ".gz"):
        if _o.path.exists(q):
            try:
                for l in (_g.open(q, "rt") if q.endswith(".gz") else open(q)):
                    yield l
            except EOFError:
                pass


ap = argparse.ArgumentParser()
ap.add_argument("run"); ap.add_argument("names", nargs="+")
ap.add_argument("--window", type=int, default=20000); ap.add_argument("--json")
a = ap.parse_args()
founders = dict(n.split("=", 1) for n in a.names)

def ropen(p):
    return open(p) if os.path.exists(p) else gzip.open(p + ".gz", "rt")

seq, tpl = {}, {}
for l in species_lines(f"{a.run}/species.tsv"):
    p = l.rstrip("\n").split("\t")
    if p[0].isdigit():
        seq[int(p[0])] = p[-1]; tpl[int(p[0])] = int(p[3])
root_of_seq = {v: k for k, v in founders.items()}
memo = {}
def lineage(s):
    path = []
    while s not in memo:
        if seq.get(s) in root_of_seq and tpl.get(s, -1) == -1:
            memo[s] = root_of_seq[seq[s]]; break
        path.append(s)
        p = tpl.get(s, -1)
        if p == -1 or p == s or len(path) > 100000:
            memo[s] = "unknown"; break
        s = p
    lab = memo[s]
    for x in path:
        memo[x] = lab
    return lab

last = None; pop = collections.Counter()
for l in ropen(f"{a.run}/popdy001.dat"):
    t, s, c = l.split(","); t = int(t)
    if t != last:
        last = t; pop = collections.Counter()
    pop[int(s)] += int(c)
n = sum(pop.values())
alive = n > 0 and "EXTINCT" not in open(f"{a.run}/log.txt").read()

rows = [l.rstrip("\n").split("\t") for l in ropen(f"{a.run}/copies.tsv")]
end = int(rows[-1][0]) if rows else 0
cp = collections.Counter()           # (maker, template) -> placed exact copies
for r in rows:
    if int(r[0]) > end - a.window and r[2] == r[3] and len(r) > 6:
        cp[(int(r[1]), int(r[2]))] += int(r[6])
tot = sum(cp.values())

out = {"run": a.run, "step": last, "alive": alive, "n": n}
lp = collections.Counter()
for s, c in pop.items():
    lp[lineage(s)] += c
out["pop"] = {k: v / n for k, v in lp.items()} if n else {}
cls = collections.Counter()
for (m, t), v in cp.items():
    key = "self" if m == t else f"{lineage(m)}->{lineage(t)}"
    cls[key] += v
out["copies"] = {k: v / tot for k, v in cls.most_common()} if tot else {}
ab = {"A->B": collections.Counter(), "B->A": collections.Counter()}
for (m, t), v in cp.items():
    lm, lt = lineage(m), lineage(t)
    if (lm, lt) == ("A", "B"): ab["A->B"][(m, t)] += v
    if (lm, lt) == ("B", "A"): ab["B->A"][(m, t)] += v
out["cross"] = {}
for k, c in ab.items():
    s = sum(c.values())
    top = c.most_common(1)[0][1] / s if s else 0
    eff = (s * s / sum(v * v for v in c.values())) if s else 0  # effective number of species pairs
    out["cross"][k] = {"copies": s, "share": s / tot if tot else 0, "pairs": len(c), "eff_pairs": round(eff, 1), "top_pair_share": round(top, 2)}
roles = []
made = collections.Counter(); bys = collections.Counter(); byl = collections.Counter(); byo = collections.Counter()
for (m, t), v in cp.items():
    made[m] += v
    if m == t: bys[t] += v
    elif lineage(m) == lineage(t): byl[t] += v
    else: byo[t] += v
for s, c in pop.most_common():
    if c / n < 0.01:
        break
    roles.append({"lin": lineage(s), "seq": seq[s], "share": round(c / n, 3), "copier": "=" in seq[s],
                  "made": made[s], "copied_self": bys[s], "copied_lineage": byl[s], "copied_other": byo[s]})
out["roles"] = roles
if a.json:
    json.dump(out, open(a.json, "w"))
print(f"== {a.run} step {last} n={n} {'alive' if alive else 'EXTINCT'}")
print("  pop    " + " ".join(f"{k}:{v:.2f}" for k, v in sorted(out['pop'].items(), key=lambda kv: -kv[1])))
print("  copies " + " ".join(f"{k}:{v:.2f}" for k, v in list(out['copies'].items())[:8]))
for k, v in out["cross"].items():
    print(f"  {k}: {v['copies']} placed ({v['share']:.2f} of all), {v['pairs']} species pairs, effective {v['eff_pairs']}, top pair {v['top_pair_share']}")
for r in roles[:12]:
    print(f"   {r['share']:.3f} {r['lin']:<7} {'=' if r['copier'] else ' '} {r['seq'][:28]:<28} made {r['made']:>6} copied: self {r['copied_self']:>5} lin {r['copied_lineage']:>5} other {r['copied_other']:>5}")
