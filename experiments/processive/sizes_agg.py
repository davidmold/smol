#!/usr/bin/env python3
"""Score the introduction-size runs (run on aida after sizes.py).

    sizes_agg.py [--scan]   (--scan first runs transition_scan.py on every finished run)

Per run: transition_scan.py --every 20000 --min 2 (species with >= 2 molecules probed; shared
cache), saved as scan.json in the run dir. Per recipient x sequence x N: runs, established
(processive share >= 10% at the last sample), persisted (any processive species at the last
sample), extinct, median final share, median step at which the share first reached 10%.
"""
import glob, json, os, statistics as st, subprocess, sys

H = os.path.expanduser("~/smol"); R = f"{H}/runs/sizes"
CACHE = "/tmp/c_sizes.json"

def finished(d):
    t = open(f"{d}/log.txt", errors="replace").read()[-4000:]
    return "Step 2700000 done" in t or "EXTINCT" in t

runs = sorted(os.path.dirname(p) for p in glob.glob(f"{R}/*/*/n*/rep*/log.txt"))
if "--scan" in sys.argv:
    if not os.path.exists(CACHE) and os.path.exists("/tmp/cs.json"):
        subprocess.run(["cp", "/tmp/cs.json", CACHE])
    for d in runs:
        if finished(d) and not os.path.exists(f"{d}/scan.json"):
            subprocess.run([sys.executable, "-I", f"{H}/tools/transition_scan.py", d, "--every", "20000",
                            "--start", "2500000", "--min", "2", "--cache", CACHE, "--json", f"{d}/scan.json"],
                           stdout=subprocess.DEVNULL, check=True)

groups = {}
for d in runs:
    if not os.path.exists(f"{d}/scan.json"):
        continue
    comm, seq, n, rep = d.split("/")[-4:]
    r = json.load(open(f"{d}/scan.json"))["rows"]
    ext = "EXTINCT" in open(f"{d}/log.txt", errors="replace").read()[-4000:]
    last = r[-1] if r else {"proc_share": 0, "proc_species": 0, "step": 0}
    t10 = next((x["step"] - 2500000 for x in r if x["proc_share"] >= 0.10), None)
    groups.setdefault((comm, seq, int(n[1:])), []).append((last["proc_share"], last["proc_species"], ext, t10))

print(f"{'recipient':10} {'seq':5} {'N':>3} {'runs':>4} {'estab':>5} {'persist':>7} {'extinct':>7} {'med_share':>9} {'med_t10':>8}")
for (c, s, n), v in sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2])):
    est = sum(1 for x in v if x[0] >= 0.10); per = sum(1 for x in v if x[1] > 0); ex = sum(1 for x in v if x[2])
    t = [x[3] for x in v if x[3] is not None]
    print(f"{c:10} {s:5} {n:>3} {len(v):>4} {est:>5} {per:>7} {ex:>7} {st.median(x[0] for x in v):>9.3f} "
          f"{(str(int(st.median(t) // 1000)) + 'k') if t else '-':>8}")
