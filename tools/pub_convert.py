#!/usr/bin/env python3
"""Convert the published spatial Stringmol run snapshots (Hickinbotham, Stepney & Hogeweg 2021,
biopaperdata/raw/runN/out1_STEP.conf, one full state every 20k steps) into our formats, so the
same analysis tools run on them.

    pub_convert.py RAW_RUN_DIR OUT_DIR [EVERY]

Writes OUT_DIR/popdy001.dat (step,species,count at each snapshot), OUT_DIR/species.tsv
(id, first step seen, -1, -1, sequence; ids assigned here, in order of first appearance) and
OUT_DIR/log.txt (one line per snapshot, so "Step N done" checks work).
"""
import os, re, sys

raw, out = sys.argv[1], sys.argv[2]
every = int(sys.argv[3]) if len(sys.argv) > 3 else 1  # keep snapshots at multiples of this
os.makedirs(out, exist_ok=True)
snaps = []
for f in os.listdir(raw):
    m = re.fullmatch(r"out1_(\d+)\.conf", f)
    if m and int(m.group(1)) % every == 0:
        snaps.append((int(m.group(1)), f))
snaps.sort()
ids = {}
pd = open(os.path.join(out, "popdy001.dat"), "w")
sp = open(os.path.join(out, "species.tsv"), "w")
lg = open(os.path.join(out, "log.txt"), "w")
for step, f in snaps:
    counts = {}
    with open(os.path.join(raw, f), errors="replace") as fh:
        for line in fh:
            # free molecules are AGENT lines; each bound pair is a REACTION block whose two
            # '###...  ID SEQ' lines give the active's and passive's species at binding
            if line.startswith("AGENT") or line.startswith("###"):
                p = line.split()
                if len(p) >= 2 and not line.startswith("###") or len(p) >= 3:
                    s = p[1] if line.startswith("AGENT") else p[2]
                    counts[s] = counts.get(s, 0) + 1
    rows = []
    for seq, c in counts.items():
        if seq not in ids:
            ids[seq] = len(ids) + 1
            sp.write(f"{ids[seq]}\t{step}\t-1\t-1\t{seq}\n")
        rows.append((ids[seq], c))
    for i, c in sorted(rows):
        pd.write(f"{step},{i},{c}\n")
    lg.write(f"Step {step} done, number of molecules is {sum(counts.values())}\n")
print(f"{out}: {len(snaps)} snapshots, {snaps[0][0]}..{snaps[-1][0]}, {len(ids)} species")
