#!/usr/bin/env python3
"""Convert a C++ Stringmol snapshot (out1_STEP.conf, e.g. from the published 2021 runs) into a
smol checkpoint, so smol can carry the run on.

    pub_to_ckpt.py OUT1_CONF OUT_DIR [--nsteps N] [--set KEY VALUE ...] [--new-rules]

Writes OUT_DIR/checkpoint_STEP.txt.gz and OUT_DIR/species.tsv. Then:

    cd OUT_DIR && smol resume checkpoint_STEP.txt.gz --reseed S

--reseed is required: smol's random generator is not the C++ one, so the C++ RNG state can't be
carried over and the checkpoint's RNG line is a placeholder. The world state carries over
exactly: every free molecule, and every bound pair mid-reaction with its strings (including
any half-written copy) and all four pointers on both strings (the C++ `irwf` fields map
one-to-one onto smol's Pointers: tape 1 = own string, tape 0 = partner).

By default the config uses the C++ partner, turn and copy-deletion rules (COMPAT_PARTNER 1,
ENCOUNTER_TURN 1, COMPAT_DELETION 1), as the published runs (v0.2.3.4) did; --new-rules leaves
them off. Parameters taken from the
snapshot: GRIDX, GRIDY, DECAY, MUTATE. Species ids are new (smol's), assigned in order of first
appearance; the C++ id of each is kept in species.tsv's third column.
"""
import argparse, gzip, os, re, sys

ap = argparse.ArgumentParser()
ap.add_argument("conf"); ap.add_argument("out")
ap.add_argument("--nsteps", type=int, default=0, help="total steps to run to (default: +1M)")
ap.add_argument("--set", nargs=2, action="append", default=[], metavar=("KEY", "VALUE"))
ap.add_argument("--new-rules", action="store_true")
a = ap.parse_args()

params = {}
# The snapshot lists every molecule in the C++ agent list's order: a free one as AGENT + GRIDPOS,
# an active one as REACTION / "### spp" / REACTION_ACTIVE (with pass_index, the position of its
# partner in this same list), a passive one as PASSIVE / REACTION_PASSIVE / "### spp" wherever it
# falls in the list.
agents = []   # dicts: kind F/A/P, seq, x, y, (A) sa, irwf, pidx, (P) sp
lines = open(a.conf, errors="replace").read().split("\n")
def grid(l):
    m = re.search(r"grid:\s+(\d+)\s+(\d+)", l)
    return int(m.group(1)), int(m.group(2))
last_hash = None
for i, l in enumerate(lines):
    p = l.split()
    if not p or l.startswith("%"):
        continue
    k = p[0]
    if k in ("EXTIT", "GRIDX", "GRIDY", "DECAY", "MUTATE", "RANDSEED", "NUMAGENTS"):
        params[k] = p[1]
    elif k == "AGENT":
        n = int(p[2])
        for j in range(n):
            q = lines[i + 1 + j].split()
            assert q[0] == "GRIDPOS", lines[i + 1 + j]
            agents.append(dict(kind="F", seq=p[1], x=int(q[1]), y=int(q[2])))
    elif k.startswith("###"):
        if agents and agents[-1]["kind"] == "P" and "sp" not in agents[-1]:
            agents[-1]["sp"] = p[2]          # closes a passive record
        else:
            last_hash = p[2]                 # opens an active record
    elif k == "REACTION_ACTIVE":
        nums = [int(v) for v in l.split("irwf:")[1].split("grid:")[0].split()]
        pidx = int(re.search(r"pass_index:\s*(-?\d+)", l).group(1))
        x, y = grid(l)
        agents.append(dict(kind="A", seq=p[2], x=x, y=y, sa=last_hash, irwf=nums, pidx=pidx))
        last_hash = None
    elif k == "REACTION_PASSIVE":
        x, y = grid(l)
        agents.append(dict(kind="P", seq=p[2], x=x, y=y))

free = [g for g in agents if g["kind"] == "F"]
reactions = []
for g in agents:
    if g["kind"] == "A":
        q = agents[g["pidx"]]
        assert q["kind"] == "P" and "sp" in q, (g, q)
        reactions.append(dict(sa=g["sa"], sp=q["sp"], astr=g["seq"], pstr=q["seq"], irwf=g["irwf"],
                              ag=(g["x"], g["y"]), pg=(q["x"], q["y"])))
if sum(g["kind"] == "P" for g in agents) != len(reactions):
    sys.exit("passive molecules without an active partner")
free = [(g["seq"], g["x"], g["y"]) for g in free]

step = int(params["EXTIT"])
gx, gy = int(params["GRIDX"]), int(params["GRIDY"])
nmol = len(free) + 2 * len(reactions)
if "NUMAGENTS" in params and int(params["NUMAGENTS"]) != nmol:
    sys.exit(f"read {nmol} molecules but NUMAGENTS is {params['NUMAGENTS']}")

ids = {}
def sid(seq):
    if seq not in ids:
        ids[seq] = len(ids) + 1
    return ids[seq]

cfg = [f"GRIDX {gx}", f"GRIDY {gy}", f"DECAY {params['DECAY']}", f"MUTATE {params['MUTATE']}",
       f"RANDSEED {params.get('RANDSEED', 1)}", f"NSTEPS {a.nsteps or step + 1_000_000}"]
if not a.new_rules:
    cfg += ["COMPAT_PARTNER 1", "ENCOUNTER_TURN 1", "COMPAT_DELETION 1"]
cfg += [f"{k} {v}" for k, v in a.set]

mols, occupied, clamped = [], set(), 0
def place(x, y):
    if not (0 <= x < gx and 0 <= y < gy):
        sys.exit(f"cell {x},{y} outside {gx}x{gy}")
    if (x, y) in occupied:
        sys.exit(f"two molecules in cell {x},{y}")
    occupied.add((x, y))
for seq, x, y in free:
    place(x, y)
    mols.append(f"{len(mols)} {x} {y} {sid(seq)} {seq} F")
for r in reactions:
    ai, pi = len(mols), len(mols) + 1
    place(*r["ag"]); place(*r["pg"])
    n = r["irwf"]
    pos, side = [], []
    for k in range(4):
        t, on_partner, on_own = n[3 * k:3 * k + 3]
        if on_partner < 0 or on_own < 0:
            clamped += 1   # the unused side of a pointer can hold a stale offset; never read
        pos += [max(on_partner, 0), max(on_own, 0)]
        side.append(t)
    a_id, p_id = sid(r["sa"]), sid(r["sp"])
    mols.append(f"{ai} {r['ag'][0]} {r['ag'][1]} {a_id} {r['astr']} A {pi} {a_id} {p_id} "
                + " ".join(map(str, pos)) + " " + " ".join(map(str, side)) + " 0e0 0")
    mols.append(f"{pi} {r['pg'][0]} {r['pg'][1]} {p_id} {r['pstr']} P {ai}")

os.makedirs(a.out, exist_ok=True)
name = os.path.join(a.out, f"checkpoint_{step}.txt.gz")
with gzip.open(name, "wt") as z:
    z.write("SMOLCKPT 1\n")
    z.write(f"STEP {step}\n")
    z.write("RNG 1 2 3 4\n")  # placeholder: resume with --reseed
    z.write(f"NEXTSPECIES {len(ids) + 1}\n")
    z.write(f"CONFIG {len(cfg)}\n" + "".join(c + "\n" for c in cfg))
    z.write(f"SLOTS {len(mols)}\nFREE 0\nMOLS {len(mols)}\n")
    z.write("".join(m + "\n" for m in mols))
with open(os.path.join(a.out, "species.tsv"), "w") as f:
    for seq, i_ in ids.items():
        f.write(f"{i_}\t{step}\t-1\t-1\t{seq}\n")
print(f"{name}: step {step}, {len(free)} free + {len(reactions)} bound pairs = {nmol} molecules, "
      f"{len(ids)} species, {gx}x{gy}" + (f", {clamped} stale pointer offsets clamped" if clamped else ""))
