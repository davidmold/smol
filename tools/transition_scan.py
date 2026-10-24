#!/usr/bin/env python3
"""Scan a run for processive replication (the classic world's second transition).

    transition_scan.py RUN_DIR [--every 100000] [--min 5] [--cap 4000] [--strand]
                       [--start 0] [--max 10**9] [--cache FILE]

Definition (24 Oct):
  processive species  in a probe self-reaction (smol-react --self, mutation off, 5000 steps)
                      it binds itself (P > 0) and makes >= 3 exact copies of itself in that
                      one reaction. Merely reaching the step limit does not count.
  detection           first sampled step at which a processive species has >= 20 molecules
  establishment       first sampled step from which processive species hold >= 10% of
                      molecules for at least 3 consecutive samples
Every species with at least --min molecules at a sampled step is probed (up to --cap per step,
most abundant first), not only the top few. Also reports cumulative reproductive activity up to
each sample: placed copies if the run logged them (OFFSPRINGLOG), else attempts.
"""
import argparse, gzip, json, os, re, subprocess, tempfile

ap = argparse.ArgumentParser(); ap.add_argument("run")
ap.add_argument("--every", type=int, default=100000); ap.add_argument("--min", type=int, default=5)
ap.add_argument("--cap", type=int, default=4000); ap.add_argument("--strand", action="store_true")
ap.add_argument("--start", type=int, default=0); ap.add_argument("--max", type=int, default=10**12)
ap.add_argument("--cache"); ap.add_argument("--json")
ap.add_argument("--react", default=os.path.expanduser("~/smol/target/release/smol-react"))
a = ap.parse_args()

def ropen(p):
    if os.path.exists(p):
        return open(p)
    return gzip.open(p + ".gz", "rt")

def lines(p):  # tolerate a gzip file cut off by a killed run
    try:
        for l in ropen(p):
            yield l
    except EOFError:
        return

pops = {}
for l in lines(f"{a.run}/popdy001.dat"):
    t = int(l[:l.index(",")])
    if t > a.max:
        break
    if t >= a.start and t % a.every == 0 and t > 0:
        _, s, c = l.split(","); pops.setdefault(t, {})[int(s)] = int(c)
need = set()
for t, p in pops.items():
    need |= {s for s, c in sorted(p.items(), key=lambda kv: -kv[1])[:a.cap] if c >= a.min}
seq = {}
for l in lines(f"{a.run}/species.tsv"):
    i = l.find("\t")
    if i > 0 and l[:i].isdigit() and int(l[:i]) in need:
        seq[int(l[:i])] = l.rstrip("\n").split("\t")[-1]

cache = json.load(open(a.cache)) if a.cache and os.path.exists(a.cache) else {}
todo = sorted({seq[s] for s in need if s in seq} - set(cache))
pat = re.compile(r"(\S+)\s+initiates with (\S+)\s+P\(bind\)=([\d.]+).*steps=(\d+)")
for i in range(0, len(todo), 400):
    chunk = todo[i:i + 400]
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        for k, sq in enumerate(chunk):
            f.write(f"R m{k} {sq}\n")
        fn = f.name
    out = subprocess.run([a.react, fn, "5000", "--self"] + (["--strand"] if a.strand else []),
                         capture_output=True, text=True).stdout.split("\n")
    os.unlink(fn)
    cur = None; res = {}
    for l in out:
        m = pat.match(l)
        if m:
            cur = m.group(1); res[cur] = [float(m.group(3)), int(m.group(4)), 0]
        elif "-> product" in l and cur and l.split()[2] == cur:
            res[cur][2] += 1
    for k, sq in enumerate(chunk):
        p, st, cp = res.get(f"m{k}", [0, 0, 0])
        cache[sq] = [p, st, cp]
if a.cache:
    json.dump(cache, open(a.cache, "w"))
proc = lambda sq: sq in cache and cache[sq][0] > 0 and cache[sq][2] >= 3

# cumulative reproductive activity
cum = {}; tot = 0; kind = None; steps = sorted(pops); j = 0
for l in lines(f"{a.run}/copies.tsv"):
    r = l.rstrip("\n").split("\t")
    t = int(r[0])
    while j < len(steps) and t > steps[j]:
        cum[steps[j]] = tot; j += 1
    if kind is None:
        kind = "placed" if len(r) > 6 else "attempted"
    tot += int(r[6]) if len(r) > 6 else int(r[4])
for s in steps[j:]:
    cum[s] = tot

print(f"== {a.run}  (cumulative {kind or '-'} copies; species with >= {a.min} molecules probed)")
print(f"{'step':>9} {'N':>7} {'probed':>6} {'cover':>5} {'proc_sp':>7} {'proc_share':>10} {'cum_copies':>12}  top processive")
det = est = None; run_len = 0; rows = []
for t in steps:
    p = pops[t]; n = sum(p.values())
    probed = [(s, c) for s, c in sorted(p.items(), key=lambda kv: -kv[1])[:a.cap] if c >= a.min and s in seq]
    pr = [(s, c) for s, c in probed if proc(seq[s])]
    share = sum(c for _, c in pr) / n if n else 0
    cover = sum(c for _, c in probed) / n if n else 0
    top = max(pr, key=lambda x: x[1]) if pr else None
    if det is None and any(c >= 20 for _, c in pr):
        det = t
    run_len = run_len + 1 if share >= 0.10 else 0
    if est is None and run_len >= 3:
        est = steps[steps.index(t) - 2]
    rows.append({"step": t, "n": n, "proc_share": share, "proc_species": len(pr), "cum": cum.get(t, 0)})
    tp = f"{top[1]:>5} {seq[top[0]][:50]}" if top else ""
    print(f"{t:>9} {n:>7} {len(probed):>6} {cover:>5.2f} {len(pr):>7} {share:>10.3f} {cum.get(t, 0):>12}  {tp}")
print(f"detection {det}  establishment {est}")
if a.json:
    json.dump({"run": a.run, "detection": det, "establishment": est, "rows": rows, "kind": kind}, open(a.json, "w"))
