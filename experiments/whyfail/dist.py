#!/usr/bin/env python3
"""Distributions behind compare.py's medians: established / large (peak >= 500, not established) /
failed lineages, by starting state. Usage: dist.py compare_rows.json"""
import json, statistics as st, sys

R = json.load(open(sys.argv[1]))
for r in R:
    c = r["copying"]
    r["state"] = "10.0M" if "s21_10.0M" in r["run"] else "9.9M"
    r["used_member"] = c.get("used_member", 0) / c["used"] if c.get("used") else None
    r["kin_share"] = c.get("kin", 0) / max(1, sum(c.get(q, 0) for q in ("kin", "used", "uses")))
    r["growth2k"] = (r["n2k"] or 0) / max(1, r["n0"] or 1)
    r["group"] = "established" if r["est"] else ("large" if r["peak"] >= 500 else "failed")
feats = ["used_member", "kin_share", "len", "selfcopies", "growth2k", "empty", "member_nb", "resident_nb", "cluster"]
for state in ("9.9M", "10.0M"):
    G = [r for r in R if r["state"] == state]
    print(f"== {state}: {len(G)} lineages with a foothold (>= 5 molecules)")
    for g in ("established", "large", "failed"):
        H = [r for r in G if r["group"] == g]
        if not H:
            continue
        print(f"  {g:11} n={len(H):3} " + "  ".join(
            f"{f}={st.median([r[f] for r in H if r[f] is not None]):.2f}" if [r for r in H if r[f] is not None] else f"{f}=-"
            for f in feats))
    F = [r for r in G if r["group"] == "failed" and r["used_member"] is not None]
    if F:
        q = sorted(r["used_member"] for r in F)
        print(f"  failed used_member: n={len(q)}  >0: {sum(x > 0 for x in q)}  >=0.4: {sum(x >= 0.4 for x in q)}  max {q[-1]:.2f}")
    print("  established/large individually:")
    for r in G:
        if r["group"] != "failed":
            print(f"    {r['group']:11} {r['run']:28} len {r['len']:3} selfcopies {r['selfcopies']:3} used_member "
                  f"{r['used_member'] if r['used_member'] is None else round(r['used_member'], 2)} kin_share {r['kin_share']:.2f} "
                  f"n0/2k {r['n0']}/{r['n2k']} peak {r['peak']} used {r['copying'].get('used', 0)}")
