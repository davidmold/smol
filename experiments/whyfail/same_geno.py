#!/usr/bin/env python3
"""Genotype or surroundings? (run on aida in ~/published/whyfail)

    same_geno.py compare_rows.json LINEAGE_JSON...

1. Same genotype, different fate: for each established or large lineage (peak >= 500) in the 9.9M
   state, find failed lineages (>= 5 molecules) elsewhere whose founder has the same sequence, or
   whose founder is one of its member sequences.
2. Who copies newcomers faithfully: in the window t5..t5+5000, the resident species that ran their
   code on a lineage member (member passive); for each, attempts and the share whose product was a
   lineage member; residents classed by whether they carry the lineage's copy core
   (the founder's last 20 symbols) - i.e. are relatives from the same resident family.
"""
import gzip, json, os, sys
from collections import Counter, defaultdict

rows = json.load(open(sys.argv[1]))
L = {}
for jf in sys.argv[2:]:
    J = json.load(open(jf))
    for x in J["lineages"]:
        L[(J["run"].split("/runs/")[-1], x["root"])] = x
R = {(r["run"], r["root"]): r for r in rows}
good = [k for k, r in R.items() if "s21_10.0M" not in k[0] and (r["est"] or r["peak"] >= 500)]
bad = [k for k, r in R.items() if "s21_10.0M" not in k[0] and not (r["est"] or r["peak"] >= 500)]
print("1. same founder genotype, different fate")
for g in good:
    x = L[g]; mseqs = {m["seq"] for m in x["members"]}
    same = [b for b in bad if L[b]["seq"] == x["seq"]]
    memb = [b for b in bad if L[b]["seq"] in mseqs and L[b]["seq"] != x["seq"]]
    print(f"  {g[0]:28} founder len {x['len']} peak {x['peak']}: failed elsewhere with the identical founder: {len(same)} "
          f"{[s[0] for s in same][:5]}; failed lineages founded by one of its top-20 member sequences: {len(memb)}")
seqcount = Counter(L[k]["seq"] for k in good + bad)
print(f"  founder sequences seen in more than one lineage: {sum(1 for s, c in seqcount.items() if c > 1)} "
      f"of {len(seqcount)} distinct founders")

def lines(p):
    for q in (p, p + ".gz"):
        if os.path.exists(q):
            try:
                for l in (gzip.open(q, "rt") if q.endswith(".gz") else open(q)):
                    yield l
            except EOFError:
                pass

print("\n2. residents that ran their code on lineage members (t5..t5+5000)")
H = os.path.expanduser("~/smol/runs")
for grp, keys in (("established/large", good), ("failed (largest 15 by attempts)", None)):
    if keys is None:
        keys = sorted(bad, key=lambda k: -R[k]["copying"].get("used", 0))[:15]
    tot = Counter()
    for k in keys:
        x = L[k]; r = R[k]; M = set(x["member_ids"]); t5 = r["t5"]; core = x["seq"][-20:]
        by = defaultdict(lambda: [0, 0])
        for l in lines(f"{H}/{k[0]}/copies.tsv"):
            p = l.split("\t"); t = int(p[0])
            if t < t5:
                continue
            if t > t5 + 5000:
                break
            a, pas, pr, att = int(p[1]), int(p[2]), int(p[3]), int(p[4])
            if pas in M and a not in M:
                by[a][0] += att; by[a][1] += att if pr in M else 0
        ids = set(by)
        seq = {}
        for l in lines(f"{H}/{k[0]}/species.tsv"):
            i = l.find("\t")
            if i > 0 and l[:i].isdigit() and int(l[:i]) in ids:
                seq[int(l[:i])] = l.rstrip("\n").rsplit("\t", 1)[-1]
        for a, (att, mem) in by.items():
            rel = "relative" if core in seq.get(a, "") else "other"
            tot[(rel, "attempts")] += att; tot[(rel, "member_products")] += mem
        top = sorted(by.items(), key=lambda kv: -kv[1][0])[:2]
        print(f"  {grp[:5]} {k[0]:28} len {x['len']:3} used {sum(v[0] for v in by.values()):4} -> member "
              f"{sum(v[1] for v in by.values()):4}; top copiers: " +
              "; ".join(f"{v[0]}x ({v[1]} member) {'REL' if core in seq.get(a, '') else 'oth'} len {len(seq.get(a, ''))} "
                        f"{seq.get(a, '?')[:40]}" for a, v in top))
    for rel in ("relative", "other"):
        a_ = tot[(rel, "attempts")]; m_ = tot[(rel, "member_products")]
        print(f"   {grp}: by {rel:8}: {a_} attempts, {m_} gave a member ({(m_ / a_ if a_ else 0):.2f})")
