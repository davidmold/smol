#!/usr/bin/env python3
"""gmt_setup.py: the tail-check gene-and-machine runs. Each gm/sN/run.conf is copied with the seed
machine and gene replaced by the tail-check machine MT (an extra `?OMSJ}` test: the partner must
carry BZFW, the gene's tail tag, or the reaction ends) and its gene GT (seed head + encode(MT)),
plus QSEARCH 1 so the test searches the whole gene. Writes gmt/sN/run.conf and gmt/jobs.txt."""
import os, subprocess
H = os.path.expanduser("~"); B = f"{H}/smol-src/target/release"; R = f"{H}/smol/runs"
M0 = "WGEWLHHRJXQQ^B?WGEW}$=?>%^F$LUBO>B^C$=?>%}"
G0 = "JTRJYUUEWKQQYHOBAOAFAIAOACADADCVALBIAJAJBGAXBQAOAFAIAOFWCZAABQAPBZBGCPCZACBEAXABAPAXBGAMCZAABQAPBZFW"
MT = "WGEWLHHRJXQQ^B?WGEW}?OMSJ}$=?>%^F$LUBO>B^C$=?>%}"
GT = G0[:16] + subprocess.run([f"{B}/smol-codon", "encode", MT], capture_output=True, text=True).stdout.strip()
jobs = open(f"{R}/gmt/jobs.txt", "w")
for s in range(301, 311):
    d = f"{R}/gmt/s{s}"; os.makedirs(d, exist_ok=True)
    out = []
    for l in open(f"{R}/gm/s{s}/run.conf"):
        if l.startswith("AGENT"):
            l = l.replace(M0, MT).replace(G0, GT)
        out.append(l)
    out.append("QSEARCH 1\n")
    open(f"{d}/run.conf", "w").write("".join(out))
    for f in ("log.txt",):
        if os.path.exists(f"{d}/{f}"): os.remove(f"{d}/{f}")
    jobs.write(f"cd {d} && {B}/smol run run.conf > log.txt 2>&1\n")
jobs.close()
print("MT", MT); print("GT", GT)
