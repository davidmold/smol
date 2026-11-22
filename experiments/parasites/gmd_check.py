#!/usr/bin/env python3
"""gmd_check.py [TS_DIR] -> pre-run check of the dummy-check control machine MD.

MD is the tail-check machine MT with its second test `?OMSJ}` (the gene must carry the tail tag
BZFW) replaced by a second copy of the head test `?WGEW}` (the gene must carry JTRJ, which the
first test already requires). MD has the same length, symbol set and number of tests as MT but
excludes no more genes than the seed machine. It separates "a recognition step that excludes
parasites" from "a longer, slower machine with an extra test".

For the seed machine M0, MT and MD (each with its gene: seed head + encode(machine)), QSEARCH on:
  loop   does M + G give both G and M back (a closed loop), and the reaction's step count
  share  share of bound gene-parasites copied, for the 5 commonest gene-parasites of every
         window in TS_DIR (default ~/smol/runs/ts_gm, the section 33 runs), median over runs.
Expected if MD is a fair control: loop closes in a similar time to MT; share close to M0's.
"""
import sys, os, re, glob, subprocess, tempfile, statistics as S
B = os.path.expanduser("~/smol-src/target/release")
TS = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/smol/runs/ts_gm")
G0 = "JTRJYUUEWKQQYHOBAOAFAIAOACADADCVALBIAJAJBGAXBQAOAFAIAOFWCZAABQAPBZBGCPCZACBEAXABAPAXBGAMCZAABQAPBZFW"
MACH = {"M0": "WGEWLHHRJXQQ^B?WGEW}$=?>%^F$LUBO>B^C$=?>%}",
        "MT": "WGEWLHHRJXQQ^B?WGEW}?OMSJ}$=?>%^F$LUBO>B^C$=?>%}",
        "MD": "WGEWLHHRJXQQ^B?WGEW}?WGEW}$=?>%^F$LUBO>B^C$=?>%}"}
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gmt_mech.py")).read()
src = src.split("tag = collections.defaultdict(list)")[0].replace("TS = sys.argv[1]", f"TS = {TS!r}")
g = {"__file__": __file__}; exec(src, g); window, copied = g["window"], g["copied"]

def gene(m):
    return G0[:16] + subprocess.run([f"{B}/smol-codon", "encode", m], capture_output=True, text=True).stdout.strip()
def loop(m, gn):
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write(f"R M {m}\nR G {gn}\n"); fn = f.name
    out = subprocess.run([f"{B}/smol-react", fn, "5000", "--strand", "--translate", "1", "--stop", "--qsearch"],
                         capture_output=True, text=True).stdout
    os.unlink(fn)
    blk = re.search(r"^M\s+initiates with G.*?(?=^\S|\Z)", out, re.M | re.S)
    if not blk: return False, None, out[:300]
    named = {"G": gn, "M": m}
    prods = [named.get(a, b) for a, b in re.findall(r"-> (?:product|passive made) (\S+)\s+(\S+)", blk.group(0))]
    st = re.search(r"steps=(\d+)", blk.group(0))
    return (gn in prods and m in prods), (int(st.group(1)) if st else None), blk.group(0).strip().splitlines()[0]

windows = []
for r in sorted(os.path.basename(x) for x in glob.glob(f"{TS}/s*")):
    for t in (200000, 600000, 1000000, 1400000, 2000000):
        if os.path.exists(f"{TS}/{r}/w{t // 1000}k/classes.tsv"):
            M, ps = window(r, t)
            if ps: windows.append((r, ps))
print(f"{len(windows)} windows of gene-parasites from {TS}")
for k, m in MACH.items():
    gn = gene(m); ok, steps, head = loop(m, gn)
    per = {}
    for r, ps in windows: per.setdefault(r, []).append(copied(m, ps))
    v = [S.median([x for x in xs if x == x]) for xs in per.values() if any(x == x for x in xs)]
    print(f"{k}: {len(m)} symbols  loop closes: {ok}  steps: {steps}  share copied (median over runs) "
          f"{S.median(v):.2f}  range {min(v):.2f}-{max(v):.2f}")
    print(f"    {head}")
