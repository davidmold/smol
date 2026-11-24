#!/usr/bin/env python3
"""pw_readout.py -> does a label change change translation? (reads ~/smol/runs/pwc2/transitions.json)

"Label change" = the sequence observation (the letters after `^F$` in the commonest machine
changed). "Change in recognition" is reserved for a measured change in behaviour. For each of the
50 transitions:

1. Where the translation-start search lands (smol-search: the `$` instruction's own computation):
   for each label (old, new) on each gene background (old gene Go, new gene Gn), the target
   position (end of the best alignment of the complemented label) and the success probability.
   A label change is "silent on a gene" if old and new labels land at the same position with
   success probabilities within 0.05.
2. Whether translation gives the working machine: for each reconstructed pair (old Mo/Go, new
   Mn/Gn, forward swap Mh/Gh, reciprocal Mr/Gr), 10 reaction trials (smol-react, seeds 1-10):
   the fraction in which the machine bound to its gene copies the gene AND makes itself.
3. Swap construction: for the forward (new label onto old pair) and reciprocal (old label onto
   new pair) swaps, the method attempted and its outcome: "exact" (tag and codons found once),
   "codons" (label codons found once; tag swapped only if exact), "no unique codons" (not
   constructible this way), or "loop fails" (constructed, but the pair does not close a loop:
   a disruptive edit that the contests leave out).
Prints a per-transition table and the counts; writes pwc2/readout.json."""
import os, sys, re, json, subprocess, tempfile, collections
H = os.path.expanduser("~"); B = f"{H}/smol-src/target/release"; OUT = f"{H}/smol/runs/pwc2"
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "pw_compete.py")).read().rsplit('\n{"find"', 1)[0]
L = {"__file__": __file__, "__name__": "lib"}; exec(src, L)
comp, enc, closes = L["comp"], L["enc"], L["closes"]

def search(pairs):
    out = subprocess.run([f"{B}/smol-search"], input="".join(f"{l} {t}\n" for l, t in pairs),
                         capture_output=True, text=True).stdout.split("\n")
    return [(int(x.split()[0]), float(x.split()[1])) for x in out if x.strip()]

def works(m, g, n=10):
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write(f"R M {m}\nR G {g}\n"); fn = f.name
    ok = 0
    for seed in range(1, n + 1):
        o = subprocess.run([f"{B}/smol-react", fn, "5000", "--strand", "--translate", "1", "--stop", "--qsearch",
                            "--seed", str(seed)], capture_output=True, text=True).stdout
        blk = re.search(r"^M\s+initiates with G.*?(?=^\S|\Z)", o, re.M | re.S)
        if blk:
            named = {"G": g, "M": m}
            p = [named.get(a, b) for a, b in re.findall(r"-> (?:product|passive made) (\S+)\s+(\S+)", blk.group(0))]
            ok += g in p and m in p
    os.unlink(fn); return ok / n

def attempt(M, G, lfrom, lto):
    Ms = M.replace(f"^F${lfrom}", f"^F${lto}", 1)
    if G.count(comp(lfrom)) == 1 and G.count(enc(lfrom)) == 1:
        g = G.replace(comp(lfrom), comp(lto)).replace(enc(lfrom), enc(lto))
        return ("exact", "ok") if closes(Ms, g) else ("exact", "loop fails")
    if G.count(enc(lfrom)) == 1:
        g = G.replace(enc(lfrom), enc(lto))
        if G.count(comp(lfrom)) == 1: g = g.replace(comp(lfrom), comp(lto))
        return ("codons", "ok") if closes(Ms, g) else ("codons", "loop fails")
    return ("-", "no unique codons")

tr = json.load(open(f"{OUT}/transitions.json")); rows = []; cnt = collections.Counter()
for k, x in enumerate(tr):
    lo, ln, Go, Gn = x["old_pw"], x["new_pw"], x["Go"], x["Gn"]
    (po_o, so_o), (pn_o, sn_o), (po_n, so_n), (pn_n, sn_n) = search([(lo, Go), (ln, Go), (lo, Gn), (ln, Gn)])
    sil_o = po_o == pn_o and abs(so_o - sn_o) <= .05; sil_n = po_n == pn_n and abs(so_n - sn_n) <= .05
    w = {"old": works(x["Mo"], Go), "new": works(x["Mn"], Gn)}
    if x.get("Gh"): w["fwd"] = works(x["Mh"], x["Gh"])
    if x.get("Gr"): w["recip"] = works(x["Mr"], x["Gr"])
    fa, ra = attempt(x["Mo"], Go, lo, ln), attempt(x["Mn"], Gn, ln, lo)
    cnt[f"fwd {fa[0]}: {fa[1]}"] += 1; cnt[f"recip {ra[0]}: {ra[1]}"] += 1
    cnt["silent on old gene"] += sil_o; cnt["silent on new gene"] += sil_n; cnt["silent on both"] += sil_o and sil_n
    r = {"k": k, "set": x["set"], "run": x["run"], "old": lo, "new": ln,
         "old_gene": {"old_label": [po_o, so_o], "new_label": [pn_o, sn_o], "silent": sil_o},
         "new_gene": {"old_label": [po_n, so_n], "new_label": [pn_n, sn_n], "silent": sil_n},
         "works": w, "fwd_attempt": fa, "recip_attempt": ra, "exact_tag_old": Go.count(comp(lo)) > 0,
         "exact_tag_new": Gn.count(comp(ln)) > 0}
    rows.append(r)
    print(f"t{k:03d} {x['set']:6} {x['run']} {lo:>5}->{ln:<5} old gene: {po_o}/{so_o:.2f} vs {pn_o}/{sn_o:.2f}"
          f"{' S' if sil_o else ' *'}  new gene: {po_n}/{so_n:.2f} vs {pn_n}/{sn_n:.2f}{' S' if sil_n else ' *'}  "
          f"works " + " ".join(f"{a}={b:.1f}" for a, b in w.items()) + f"  fwd {fa[1]}  recip {ra[1]}")
n = len(tr)
print(f"\n{n} transitions. Label change silent (same start position, success prob within 0.05):")
print(f"  on the old gene: {cnt['silent on old gene']}/{n}; on the new gene: {cnt['silent on new gene']}/{n}; on both: {cnt['silent on both']}/{n}")
print(f"Exact complement of the label in its own gene: old pairs {sum(r['exact_tag_old'] for r in rows)}/{n}, new pairs {sum(r['exact_tag_new'] for r in rows)}/{n}")
for key in ("old", "new", "fwd", "recip"):
    v = [r["works"][key] for r in rows if key in r["works"]]
    if v: print(f"Working machine made, {key:5} pairs: mean {sum(v) / len(v):.2f} of trials; all 10 trials in {sum(1 for a in v if a == 1)}/{len(v)}")
print("Swap attempts:"); [print(f"  {a}: {b}") for a, b in sorted(cnt.items()) if a.startswith(("fwd", "recip"))]
json.dump(rows, open(f"{OUT}/readout.json", "w"), indent=1)
