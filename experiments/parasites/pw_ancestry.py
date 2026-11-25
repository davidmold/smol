#!/usr/bin/env python3
"""pw_ancestry.py [k ...] -> are the before and after pairs of a label transition one lineage?

species.tsv records, for every new species, the copier (pa) and the template (pp): a copied gene's
pp is the gene it was copied from, a translated machine's pp is the gene it was translated from.
So a gene's ancestry is its chain of pp links. For each transition k (pwc2/transitions.json):
- the template chains of the old gene Go and the new gene Gn back to the start;
- whether Go is an ancestor of Gn (one lineage) or, if not, their most recent common ancestor
  and when it was born (competing branches since then);
- along Gn's chain back to the common ancestor: when the new label's codons first appear (the
  label change) and how many other letters changed on that branch before and after it.
Default: all transitions with both swaps, plus any given on the command line."""
import sys, os, json, subprocess
H = os.path.expanduser("~"); R = f"{H}/smol/runs"; B = f"{H}/smol-src/target/release"
tr = json.load(open(f"{R}/pwc2/transitions.json"))
enc = lambda s: subprocess.run([f"{B}/smol-codon", "encode", s], capture_output=True, text=True).stdout.strip()
def lev(a, b):
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i] + [0] * len(b)
        for j, cb in enumerate(b, 1): cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb))
        prev = cur
    return prev[-1]
ks = [int(a) for a in sys.argv[1:]] or [k for k, x in enumerate(tr) if x.get("Gh") and x.get("Gr")]
cache = {}
for k in ks:
    x = tr[k]; d = f"{R}/{x['set']}/{x['run']}"
    if d not in cache:
        sp = {}; byseq = {}
        for l in open(f"{d}/species.tsv"):
            p = l.rstrip("\n").split("\t")
            if len(p) < 5 or not p[0].isdigit(): continue
            i = int(p[0]); sp[i] = (int(p[1]), int(p[2]), int(p[3]), p[4]); byseq.setdefault(p[4], i)
        cache = {d: (sp, byseq)}
    sp, byseq = cache[d]
    def chain(i):
        out = []
        while i in sp and i not in out:
            out.append(i); i = sp[i][2]
        return out
    go, gn = byseq.get(x["Go"]), byseq.get(x["Gn"])
    if go is None or gn is None:
        print(f"t{k:03d}: gene not found in species.tsv"); continue
    co, cn = chain(go), chain(gn); so = set(co)
    mrca = next((i for i in cn if i in so), None)
    one = go in set(cn)
    print(f"t{k:03d} {x['set']} {x['run']} {x['old_pw']}->{x['new_pw']} (windows {x['t'] // 1000}k -> {x['t2'] // 1000}k)")
    print(f"   Go species {go} born {sp[go][0]}; Gn species {gn} born {sp[gn][0]}; "
          + ("Go IS an ancestor of Gn (one lineage)" if one else
             f"competing branches: common ancestor {mrca} born {sp[mrca][0] if mrca else '?'}"))
    stop = cn.index(mrca) if mrca in cn else len(cn) - 1
    branch = list(reversed(cn[:stop + 1]))           # from the common ancestor (or Go) to Gn
    lab = enc(x["new_pw"]); first = next((i for i in branch if lab in sp[i][3]), None)
    print(f"   branch from {'Go' if one else 'common ancestor'} to Gn: {len(branch) - 1} copying steps; "
          f"total edits {lev(sp[branch[0]][3], x['Gn'])}")
    if first is not None:
        j = branch.index(first)
        print(f"   new label's codons first on species {first} (born {sp[first][0]}), step {j} of {len(branch) - 1}; "
              f"edits before it on the branch {lev(sp[branch[0]][3], sp[first][3])}, after it {lev(sp[first][3], x['Gn'])}")
    else:
        print("   new label's codons not found on the branch (codons differ; label from translation)")
