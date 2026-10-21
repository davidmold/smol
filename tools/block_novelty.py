# Per-block novelty from checkpoints, and the small runs computed the same way (one snapshot per 100k)
import sys,os,gzip,glob,collections,statistics
sys.path.insert(0,'/home/davidmold/stringmol')
from families import load_species, read_popdy, cluster
import numpy as np
from rapidfuzz.distance import Levenshtein
from rapidfuzz.process import cdist
T=0.25
def score(prev,cur,seq):
    if not prev or not cur: return None
    cur=sorted(cur);prev=sorted(prev)
    dmin=cdist([seq[s] for s in cur],[seq[s] for s in prev],scorer=Levenshtein.normalized_distance,dtype=np.float32,workers=-1).min(axis=1)
    new=[s for s,d in zip(cur,dmin) if d>T]
    nf=len(set(cluster(new,seq,T).values())) if new else 0
    return nf,len(set(cluster(cur,seq,T).values()))
def abundant(c):
    tot=sum(c.values()); lim=max(20,0.01*tot)
    return {s for s,n in c.items() if n>=lim}
D='/home/davidmold/smol/runs/big/classic_400x300'
seq={}
snaps={}
for f in sorted(glob.glob(D+'/checkpoint_*.txt.gz'),key=lambda f:int(f.split('_')[-1].split('.')[0])):
    st=int(f.split('_')[-1].split('.')[0]); B=collections.defaultdict(collections.Counter); mols=False
    for l in gzip.open(f,'rt'):
        if l.startswith('MOLS'): mols=True; continue
        if not mols: continue
        p=l.split()
        if len(p)<5: continue
        x,y,sp=int(p[1]),int(p[2]),int(p[3]); seq[sp]=p[4]
        B[(x//100,y//100)][sp]+=1
    snaps[st]=B
steps=sorted(snaps)
print('BIG, per 100x100 block, snapshot every 100k:')
for a,b in zip(steps,steps[1:]):
    r=[score(abundant(snaps[a][k]),abundant(snaps[b][k]),seq) for k in sorted(snaps[b])]
    r=[x for x in r if x]
    print(f' {a//1000}->{b//1000}k  blocks {len(r)}  families/block {statistics.mean(x[1] for x in r):.1f}  new/block {statistics.mean(x[0] for x in r):.2f}')
print('SMALL pr0 runs, same snapshot method:')
for w in (8,9):
  nn=[];ff=[]
  for d in sorted(glob.glob('/home/davidmold/smol/runs/pr0/seed*')):
    pop=read_popdy(d,1000)
    a,b=w*100000,(w+1)*100000
    if a in pop and b in pop and sum(pop[b].values())>2:
        s,_=load_species(d)
        x=score(abundant(pop[a]),abundant(pop[b]),s)
        if x: nn.append(x[0]); ff.append(x[1])
  print(f' {w*100}->{(w+1)*100}k  runs {len(nn)}  families/run {statistics.mean(ff):.1f}  new/run {statistics.mean(nn):.2f}')
print('BIG, per block, "new" only if unlike anything abundant in ANY block at the previous snapshot:')
for a,b in zip(steps,steps[1:]):
    prevall=set().union(*[abundant(snaps[a][k]) for k in snaps[a]])
    r=[score(prevall,abundant(snaps[b][k]),seq) for k in sorted(snaps[b])]
    r=[x for x in r if x]
    curall=set().union(*[abundant(snaps[b][k]) for k in snaps[b]])
    g=score(prevall,curall,seq)
    print(f' {a//1000}->{b//1000}k  new/block {statistics.mean(x[0] for x in r):.2f}   whole grid: families {g[1]} new {g[0]}')
