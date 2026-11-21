#!/usr/bin/env python3
"""gm_linkage.py RUN_DIR... -> do machines favour their own gene?

A machine species m is born by translation of a gene g_m (its parent that is not a machine: the
template, found in species.tsv as the parent whose sequence has no opcodes... in practice the
longer, mostly-letter parent). Over 1M-2M, from copies.tsv: of all exact gene copies machine m
makes, the share that are copies of g_m ("own-gene share"), against the share of g_m among all
exact copies made by all machines in the same period ("chance"). Linkage = own / chance:
1 means a machine copies its own gene no more often than any machine does, so a better machine
does not help its own gene and selection cannot act on machine changes."""
import sys, collections, statistics as S
OPS = set("$>^?=%}+")
for d in sys.argv[1:]:
    seq, par = {}, {}
    for l in open(f"{d}/species.tsv"):
        p = l.rstrip("\n").split("\t")
        if p[0].isdigit() and len(p) >= 5:
            seq[int(p[0])] = p[-1]; par[int(p[0])] = (int(p[2]), int(p[3]))
    def is_machine(s): return sum(ch in OPS for ch in seq.get(s, "")) >= 6
    own = collections.Counter(); allc = collections.Counter(); genecopies = collections.Counter()
    for l in open(f"{d}/copies.tsv"):
        p = l.split("\t")
        if len(p) < 5 or not p[0].isdigit(): continue
        st, c, t, pr, n = int(p[0]), int(p[1]), int(p[2]), int(p[3]), int(p[4])
        if st <= 1000000 or pr != t or is_machine(t) or not is_machine(c): continue
        allc[c] += n; genecopies[t] += n
        a, b = par.get(c, (-1, -1))
        gm = [x for x in (a, b) if x in seq and not is_machine(x)]
        if gm and t == gm[0]: own[c] += n
    tot = sum(genecopies.values())
    rows = []
    for c, n in allc.most_common(20):
        a, b = par.get(c, (-1, -1)); gm = [x for x in (a, b) if x in seq and not is_machine(x)]
        if not gm or n < 200: continue
        ownshare = own[c] / n; chance = genecopies[gm[0]] / tot
        rows.append((n, ownshare, chance))
    if rows:
        ratio = [o / ch for _, o, ch in rows if ch > 0]
        print(f"{d}: {len(rows)} commonest machines; own-gene share median {S.median(o for _, o, _ in rows):.3f}, "
              f"chance median {S.median(ch for _, _, ch in rows):.3f}, linkage (own/chance) median {S.median(ratio):.2f} "
              f"[{min(ratio):.2f}-{max(ratio):.2f}]")
