#!/usr/bin/env python3
"""confirm_setup.py: the confirmation round for the tail-check arms race (plan: CONFIRM_PLAN.md).

New runs, seeds 311-340, the section 33 layout (gm/s301/run.conf with RANDSEED changed), QSEARCH 1:
  cf_gmt  tail-check machine MT + its gene     (the result being confirmed)
  cf_gmq  seed machine M0 + seed gene          (same rules, no extra test)
  cf_gmd  dummy-check machine MD + its gene    (extra test of the same length that excludes
                                                nothing new: the head test twice)
Extensions to 6M of the original gmt and gmq runs (seeds 301-310), resumed from checkpoint_2000000
into gmt6/sN and gmq6/sN (the 2M originals are left untouched).
Writes runs/confirm_jobs.txt: extensions first (longest), then the new runs interleaved by seed."""
import os, subprocess
H = os.path.expanduser("~"); B = f"{H}/smol-src/target/release"; R = f"{H}/smol/runs"
T = f"{H}/smol-src/tools"
M0 = "WGEWLHHRJXQQ^B?WGEW}$=?>%^F$LUBO>B^C$=?>%}"
G0 = "JTRJYUUEWKQQYHOBAOAFAIAOACADADCVALBIAJAJBGAXBQAOAFAIAOFWCZAABQAPBZBGCPCZACBEAXABAPAXBGAMCZAABQAPBZFW"
MT = "WGEWLHHRJXQQ^B?WGEW}?OMSJ}$=?>%^F$LUBO>B^C$=?>%}"
MD = "WGEWLHHRJXQQ^B?WGEW}?WGEW}$=?>%^F$LUBO>B^C$=?>%}"
enc = lambda m: subprocess.run([f"{B}/smol-codon", "encode", m], capture_output=True, text=True).stdout.strip()
SETS = {"cf_gmt": (MT, G0[:16] + enc(MT)), "cf_gmq": (M0, G0), "cf_gmd": (MD, G0[:16] + enc(MD))}
assert SETS["cf_gmq"][1] == G0[:16] + enc(M0)
base = open(f"{R}/gm/s301/run.conf").read().splitlines(keepends=True)
jobs = []
for s in range(301, 311):
    for src, dst in (("gmt", "gmt6"), ("gmq", "gmq6")):
        d = f"{R}/{dst}/s{s}"
        if os.path.exists(f"{d}/log.txt"): continue
        jobs.append(f"python3 {T}/ckpt_edit.py {R}/{src}/s{s} checkpoint_2000000.txt.gz {d} --set NSTEPS 6000000 > /dev/null 2>&1"
                    f" && cd {d} && {B}/smol resume checkpoint_2000000.txt.gz > log.txt 2>&1")
for s in range(311, 341):
    for name, (m, g) in SETS.items():
        d = f"{R}/{name}/s{s}"; os.makedirs(d, exist_ok=True)
        out = []
        for l in base:
            if l.startswith("RANDSEED"): l = f"RANDSEED {s}\n"
            if l.startswith("AGENT"): l = l.replace(M0, m).replace(G0, g)
            out.append(l)
        out.append("QSEARCH 1\n")
        open(f"{d}/run.conf", "w").write("".join(out))
        jobs.append(f"cd {d} && {B}/smol run run.conf > log.txt 2>&1")
open(f"{R}/confirm_jobs.txt", "w").write("\n".join(jobs) + "\n")
for k, (m, g) in SETS.items(): print(k, m, g)
print(len(jobs), "jobs")
