#!/usr/bin/env python3
"""gm_scan.py [GENE MACHINE] -> why is the machine frozen? Mutational scan of a working gene.

Every single-letter substitution in the gene's coded part (from the start of the encoded machine
to the end) is made, and the machine M translates the mutant gene g' (smol-react --strand
--translate 1 --stop, mutation off). The new machine m' is:
  synonymous   m' == M
  works        m' != M, and m' bound to g' copies g' and makes m' again (a closed loop)
  broken       m' != M and the loop does not close
  no product   M makes no translation of g'
Also: for the codon table itself (smol-codon table), the share of single-letter codon changes that
keep the same symbol, and that keep the symbol's kind (letter vs opcode).
Default gene and machine: the section 8 seed. Coded part: where `smol-codon encode M` sits in G."""
import sys, os, re, subprocess, tempfile, collections
B = os.path.expanduser("~/smol-src/target/release")
M0 = "WGEWLHHRJXQQ^B?WGEW}$=?>%^F$LUBO>B^C$=?>%}"
G0 = "JTRJYUUEWKQQYHOBAOAFAIAOACADADCVALBIAJAJBGAXBQAOAFAIAOFWCZAABQAPBZBGCPCZACBEAXABAPAXBGAMCZAABQAPBZFW"
G, M = (sys.argv[1], sys.argv[2]) if len(sys.argv) > 2 else (G0, M0)
L = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

def react(m, g):
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write(f"R M {m}\nR G {g}\n"); n = f.name
    out = subprocess.run([f"{B}/smol-react", n, "5000", "--strand", "--translate", "1", "--stop"],
                         capture_output=True, text=True).stdout
    os.unlink(n)
    blk = re.search(r"^M\s+initiates with G.*?(?=^\S|\Z)", out, re.M | re.S)
    if not blk: return []
    named = {"G": g, "M": m}
    return [named.get(a, b) for a, b in re.findall(r"-> (?:product|passive made) (\S+)\s+(\S+)", blk.group(0))]

enc = subprocess.run([f"{B}/smol-codon", "encode", M], capture_output=True, text=True).stdout.strip()
start = G.find(enc) if enc and enc in G else None
if start is None:      # evolved gene: the coded part is the last 2*len(M) letters
    start = len(G) - 2 * len(M)
print(f"gene {len(G)} letters, machine {len(M)} symbols, coded part from {start}")
base = react(M, G)
print("check: M with G gives", "copy+M" if (G in base and M in base) else base)

res = collections.Counter(); bypos = collections.defaultdict(collections.Counter); works = []
for i in range(start, len(G)):
    for c in L:
        if c == G[i]: continue
        g2 = G[:i] + c + G[i + 1:]
        p = react(M, g2)
        ms = [x for x in p if x != g2]
        if not ms: k = "no product"
        elif M in ms: k = "synonymous"
        else:
            m2 = ms[0]; p2 = react(m2, g2)
            k = "works" if (g2 in p2 and m2 in p2) else "broken"
            if k == "works": works.append((i, c, m2))
        res[k] += 1; bypos[(i - start) % 2][k] += 1
tot = sum(res.values())
print("\nsingle-letter substitutions in the coded part:", tot)
for k in ("synonymous", "works", "broken", "no product"):
    print(f"  {k:12} {res[k]:5}  {res[k] / tot:.3f}")
for pos in (0, 1):
    t = sum(bypos[pos].values())
    print(f"  codon letter {pos + 1}: " + "  ".join(f"{k} {bypos[pos][k] / t:.2f}" for k in ("synonymous", "works", "broken", "no product")))
print("\nexamples of working new machines:")
for i, c, m2 in works[:8]:
    print(f"  pos {i} -> {c}: {m2}")

# codon table robustness
tab = subprocess.run([f"{B}/smol-codon", "table"], capture_output=True, text=True).stdout
code = {}
for l in tab.splitlines():
    p = l.split()
    if len(p) > 1:
        for cd in p[1:]: code[cd] = p[0]
same = kind = n = 0
for cd, sym in code.items():
    for j in (0, 1):
        for c in L:
            if c == cd[j]: continue
            cd2 = cd[:j] + c + cd[j + 1:]
            if cd2 not in code: continue
            n += 1; s2 = code[cd2]
            same += s2 == sym; kind += s2.isalpha() == sym.isalpha()
print(f"\ncodon table: {len(code)} codons, {len(set(code.values()))} symbols; single-letter codon changes that keep "
      f"the symbol {same / n:.3f}, that keep its kind (letter/opcode) {kind / n:.3f}")
