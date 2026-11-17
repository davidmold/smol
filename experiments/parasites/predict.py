#!/usr/bin/env python3
"""predict.py -> (v2: 10 trials per pair; fastest-parasite measure added) does the parasite pool's measured exploitation of A vs D predict who wins?

For every pool used in ratio.py (snapshot pools) and evolve.py (evolved pools), each parasite is
reacted with A and with D in smol-react (5000 steps, --compat). Exploitation of a host:
  E1 = sum over parasites of pool fraction x P(bind) x [host makes an exact copy]
  E2 = E1 with each term weighted by (host's own self-copy steps / steps to copy the parasite),
       i.e. how many parasite copies the host makes in the time it takes to copy itself.
  E3 = max over parasites of P(bind) x success rate x (host self-copy steps / parasite copy steps):
       the fastest-growing parasite on that host (pre-stated hypothesis: the fastest one takes over).
Each host-parasite pair and each host's self-copy is reacted in its own file over 10 seeds (--seed),
since the generator state carried between pairs occasionally changes the outcome.
Prediction: D does better than its no-parasite baseline when E(D) < E(A).
Outcome: median D share of A+D at 30k over 10 seeds, minus baseline (compete m0p0).
"""
import os, re, glob, subprocess, tempfile, statistics, collections
H = os.path.expanduser("~")
R = f"{H}/smol-src/target/release/smol-react"
PAIRS = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(f"{H}/smol/runs/compete/pairs.tsv")}

TRIALS = range(1, 11)
_cache = {}
def pair(h, p):
    """(P(bind), fraction of trials where h makes an exact copy of p, median steps of those trials)."""
    if (h, p) in _cache: return _cache[(h, p)]
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write(f"R h {h}\n" + (f"R p {p}\n" if p != h else "")); n = f.name
    pb, ok, st = 0.0, 0, []
    for t in TRIALS:
        out = subprocess.run([R, n, "5000", "--compat", "--seed", str(t)] + (["--self"] if p == h else []),
                             capture_output=True, text=True).stdout
        cur = None
        for l in out.splitlines():
            m = re.match(r"(\w+)\s+initiates with (\w+)\s+P\(bind\)=([\d.]+)\s+copier=(\w+)\s+copied=(\w+)\s+steps=(\d+)", l)
            if m:
                cur = (m.group(4), m.group(5), float(m.group(3)), int(m.group(6))); continue
            m = re.match(r"\s+-> (?:product|passive made) (\w+)", l)
            want = ("h", "h") if p == h else ("h", "p")
            if m and cur and (cur[0], cur[1]) == want and m.group(1) == want[1]:
                pb = cur[2]; ok += 1; st.append(cur[3]); break
        else:
            pass
        if not pb:
            mm = re.search(r"h\s+initiates with " + ("h" if p == h else "p") + r"\s+P\(bind\)=([\d.]+)", out)
            if mm: pb = float(mm.group(1))
    os.unlink(n)
    r = (pb, ok / len(TRIALS), statistics.median(st) if st else 0)
    _cache[(h, p)] = r; return r

def exploit(A, D, par, cnt):
    tot = sum(cnt); e = {}
    for h, H_ in (("A", A), ("D", D)):
        _, sok, sst = pair(H_, H_)
        sst = sst or 5000
        e1 = e2 = e3 = 0.0
        for p, n in zip(par, cnt):
            pb, ok, st = pair(H_, p)
            if ok:
                rate = pb * ok * sst / max(st, 1)
                e1 += n / tot * pb * ok; e2 += n / tot * rate; e3 = max(e3, rate)
        e[h] = (e1, e2, e3)
    return e

def share(d, A, D):
    seq = {}
    for l in open(f"{d}/species.tsv"):
        p = l.rstrip("\n").split("\t")
        if p[0].isdigit(): seq[int(p[0])] = p[-1]
    a = b = 0
    for l in open(f"{d}/popdy001.dat"):
        s, sp, n = (int(v) for v in l.split(","))
        if s == 30000 and sp in seq:
            if seq[sp] == A: a += n
            elif seq[sp] == D: b += n
    return b / (a + b) if a + b else None

def outcome(pattern, A, D):
    v = [x for x in (share(d, A, D) for d in glob.glob(pattern)) if x is not None]
    return statistics.median(v) if v else None

base = {run: outcome(f"{H}/smol/runs/compete/{run}/m0p0/s*", p[2], p[3]) for run, p in PAIRS.items()}
pools = []
for l in open(f"{H}/smol/runs/ratio/chosen.tsv"):
    q = l.rstrip("\n").split("\t")
    pools.append(("snap", q[0], f"src{q[1]}", q[6].split(","), [int(x) for x in q[7].split(",")],
                  f"{H}/smol/runs/ratio/{q[0]}/src{q[1]}/s*"))
for l in open(f"{H}/smol/runs/evolve/pools.tsv"):
    q = l.rstrip("\n").split("\t")
    if int(q[5]) < 200 or not q[8]: continue
    pools.append(("evo", q[0], f"on{q[1]}{q[2]}_{int(q[3]) // 1000}k", q[8].split(","), [int(x) for x in q[9].split(",")],
                  f"{H}/smol/runs/evolve/{q[0]}/comp{q[1]}{q[2]}_{int(q[3]) // 1000}k/s*"))

rows = []
print(f"{'set':5} {'run':5} {'pool':12} {'E1(A)':>6} {'E1(D)':>6} {'E2(A)':>6} {'E2(D)':>6} {'E3(A)':>6} {'E3(D)':>6}  {'Dshare':>6} {'base':>5}  pred1 pred2 pred3")
for kind, run, name, par, cnt, pat in pools:
    A, D = PAIRS[run][2], PAIRS[run][3]
    o = outcome(pat, A, D)
    if o is None: continue
    e = exploit(A, D, par, cnt)
    dwins = o > base[run]
    p1, p2, p3 = ((e["D"][k] < e["A"][k]) == dwins for k in range(3))
    rows.append((kind, run, name, e, o, dwins, p1, p2, p3))
    print(f"{kind:5} {run:5} {name:12} {e['A'][0]:6.2f} {e['D'][0]:6.2f} {e['A'][1]:6.2f} {e['D'][1]:6.2f} {e['A'][2]:6.2f} {e['D'][2]:6.2f}  {o:6.2f} {base[run]:5.2f}  "
          f"{'ok' if p1 else 'MISS':5} {'ok' if p2 else 'MISS':5} {'ok' if p3 else 'MISS'}", flush=True)
for lab, sel in (("all", lambda r: True), ("snapshot pools", lambda r: r[0] == "snap"), ("evolved pools", lambda r: r[0] == "evo")):
    rr = [r for r in rows if sel(r)]
    print(f"{lab}: E1 predicts the winner in {sum(r[6] for r in rr)}/{len(rr)}, E2 in {sum(r[7] for r in rr)}/{len(rr)}, E3 in {sum(r[8] for r in rr)}/{len(rr)}")
