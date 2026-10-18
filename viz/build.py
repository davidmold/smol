#!/usr/bin/env python3
"""Build a self-contained viewer page for a Stringmol run directory.
Usage: build.py RUN_DIR [OUT.html] [EVERY] [FROM] [TO]
(default OUT = RUN_DIR/viewer.html, EVERY = 500 steps, whole run; FROM/TO limit the steps shown)

Molecules are coloured by what they do, never by what they look like:
  1 replicator    copies others and is copied
  2 parasite      is copied, copies nobody
  3 copier only   copies others, is never copied
  4 unclassified  no copy events seen (or no data)
Roles come from RUN_DIR/roles.tsv (sequence <tab> role name) if present, else from the
simulator's copy log RUN_DIR/copies.tsv (step, active, passive, product, count): a species
counts as copying / being copied once it has done so at least MINEV times in the run.
Each frame: one byte per grid cell. 0 = empty; 1..4 = role of a non-top species;
5+ = one of the top-N lineages."""
import sys, os, re, json, base64, collections
run = sys.argv[1].rstrip('/')
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(run, 'viewer.html')
every = int(sys.argv[3]) if len(sys.argv) > 3 else 500
lo = int(sys.argv[4]) if len(sys.argv) > 4 else 0
hi = int(sys.argv[5]) if len(sys.argv) > 5 else 10**12
NTOP = 12

MINEV = 3
ROLE = {'replicator': 1, 'parasite': 2, 'copier': 3, 'inert': 4, 'unclassified': 4}
role_of = {}            # sequence -> 1..4

def load_roles():
    rp = os.path.join(run, 'roles.tsv')
    if os.path.exists(rp):
        for line in open(rp):
            a = line.rstrip('\n').split('\t')
            if len(a) >= 2: role_of[a[0]] = ROLE.get(a[1], 4)
        return 'roles.tsv'
    cp = os.path.join(run, 'copies.tsv'); sp = os.path.join(run, 'species.tsv')
    if not (os.path.exists(cp) and os.path.exists(sp)):
        return None
    seqs = {}
    for line in open(sp):
        a = line.rstrip('\n').split('\t')
        if len(a) >= 5: seqs[int(a[0])] = a[4]
    copier = collections.Counter(); copied = collections.Counter()
    for line in open(cp):
        a = line.split('\t')
        if len(a) < 5: continue
        act, pas, prod, n = int(a[1]), int(a[2]), int(a[3]), int(a[4])
        if prod == pas: copier[act] += n; copied[pas] += n      # active copied the passive
        elif prod == act: copier[act] += n; copied[act] += n    # active copied itself
    for i, q in seqs.items():
        c, d = copier[i] >= MINEV, copied[i] >= MINEV
        role_of[q] = 1 if c and d else 2 if d else 3 if c else 4
    return 'copies.tsv'

def cls(s):
    return role_of.get(s, 4)

def slim_frames():
    # slim output: species.tsv (id -> sequence) + grid.snap (zlib frames of uint32 species ids)
    import struct, zlib
    seqs = {}
    for line in open(os.path.join(run, 'species.tsv')):
        a = line.rstrip('\n').split('\t')
        if len(a) >= 5: seqs[int(a[0])] = a[4]
    out = []; gx = gy = None
    with open(os.path.join(run, 'grid.snap'), 'rb') as f:
        while True:
            h = f.read(16)
            if len(h) < 16: break
            st, gx, gy, n = struct.unpack('<4I', h); z = f.read(n)
            if len(z) < n: break                     # partial frame while the run is writing
            if st % every or not lo <= st <= hi: continue
            ids = struct.unpack('<%dI' % (gx * gy), zlib.decompress(z))
            out.append((st, [(i % gx, i // gx, seqs.get(v & 0x7fffffff, '')) for i, v in enumerate(ids) if v]))
    return out, gx, gy

src = load_roles()
print('roles from', src or 'nowhere: everything unclassified')
steps = sorted(int(m.group(1)) for f in os.listdir(run) if (m := re.match(r'out1_(\d+)\.conf$', f)))
steps = [s for s in steps if s % every == 0 and lo <= s <= hi]
frames = []
gx = gy = None
abund = collections.Counter(); first = {}
slim = os.path.exists(os.path.join(run, 'grid.snap'))
if slim:
    loaded, gx, gy = slim_frames(); steps = []
else:
    loaded = None
for s in (steps if not slim else []):
    mols = []; seq = None
    for line in open(os.path.join(run, f'out1_{s:05d}.conf' if s < 100000 else f'out1_{s}.conf')):
        if line.startswith('GRIDX'): gx = int(line.split()[1])
        elif line.startswith('GRIDY'): gy = int(line.split()[1])
        elif line.startswith('AGENT'): seq = line.split()[1]
        elif line.startswith('GRIDPOS') and seq is not None:
            _, x, y = line.split(); mols.append((int(x), int(y), seq)); seq = None
        # bound pairs are written as reaction blocks; the '#####' line carries the
        # active molecule's species sequence (its REACTION_ACTIVE string is mid-copy)
        elif line.startswith('###############'): seq = line.split()[2]
        elif line.startswith('REACTION_ACTIVE') or line.startswith('REACTION_PASSIVE'):
            p = line.split(); gi = p.index('grid:'); x, y = int(p[gi + 1]), int(p[gi + 2])
            q = seq if p[0] == 'REACTION_ACTIVE' else p[2]
            mols.append((x, y, q)); seq = None
    frames.append((s, mols))
    if not mols: break          # extinct: stop at the first empty frame
    for _, _, q in mols:
        abund[q] += 1; first.setdefault(q, s)
if slim:
    for s, mols in loaded:
        frames.append((s, mols))
        for _, _, q in mols:
            abund[q] += 1; first.setdefault(q, s)
        if not mols: break
top = [q for q, _ in abund.most_common(NTOP)]
topidx = {q: i + 5 for i, q in enumerate(top)}
enc, series = [], []
for s, mols in frames:
    g = bytearray(gx * gy)
    c = collections.Counter(); tc = collections.Counter()
    for x, y, q in mols:
        k = cls(q); c[k] += 1
        v = topidx.get(q, k); g[y * gx + x] = v
        if v >= 5: tc[v - 5] += 1
    enc.append(base64.b64encode(bytes(g)).decode())
    series.append({'t': s, 'rep': c[1], 'par': c[2], 'cop': c[3], 'unk': c[4],
                   'top': [tc[i] for i in range(len(top))]})
data = json.dumps({'roles': src, 'name': os.path.basename(os.path.abspath(run)), 'gx': gx, 'gy': gy, 'every': every,
           'species': [{'seq': q, 'cls': cls(q), 'len': len(q), 'first': first[q]} for q in top],
           'series': series, 'frames': enc}, separators=(',', ':'))
tpl = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'viewer.html')).read()
open(out, 'w').write(tpl.replace('/*DATA*/', data.replace('</', '<\\/')))
print(len(frames), 'frames', gx, 'x', gy, '->', out, os.path.getsize(out) // 1024, 'KB')
