import os, sys, collections
raw = os.path.expanduser("~/published/bio/biopaperdata/raw")
out = sys.argv[1]; K = int(sys.argv[2])
os.makedirs(out, exist_ok=True)
for r in range(1, 9):
    for step in range(100000, 2000001, 100000):
        f = f"{raw}/run{r}/out1_{step}.conf"
        if not os.path.exists(f): continue
        c = collections.Counter()
        for l in open(f, errors="replace"):
            if l.startswith("AGENT"):
                p = l.split(); c[p[1]] += 1
            elif l.startswith("###"):
                p = l.split()
                if len(p) >= 3: c[p[2]] += 1
        with open(f"{out}/r{r}_{step}.txt", "w") as o:
            for i, (s, n) in enumerate(c.most_common(K)):
                o.write(f"R s{i+1} {s}\n")
