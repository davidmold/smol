import sys,os,collections
D=sys.argv[1]; W=int(sys.argv[2]) if len(sys.argv)>2 else 100000
seqs={}
for l in open(D+'/species.tsv'):
    p=l.rstrip('\n').split('\t')
    if p[0].isdigit(): seqs[int(p[0])]=p[-1]
f=open(D+'/copies.tsv','rb'); f.seek(0,2); sz=f.tell(); f.seek(max(0,sz-80_000_000)); f.readline()
rows=[]
for l in f:
    p=l.decode().split('\t')
    if len(p)<5 or not p[0].isdigit(): continue
    rows.append((int(p[0]),int(p[1]),int(p[2]),int(p[3]),int(p[4])))
last=rows[-1][0]; rows=[r for r in rows if r[0]>last-W]
cop=lambda i:'=' in seqs.get(i,'')
tot=collections.Counter(); pair=collections.Counter(); ex=0; allev=0
for s,a,pa,pr,n in rows:
    allev+=n
    if pr!=pa: continue
    ex+=n
    k=('copier' if cop(a) else 'NOcode')+'->'+('copier' if cop(pa) else 'passenger')+('(self)' if a==pa else '')
    tot[k]+=n; pair[(a,pa)]+=n
print('==',D.split('/')[-1],'steps',last-W,'-',last,' events',allev,' exact copies',ex,'(%.0f%%)'%(100*ex/allev))
for k,v in tot.most_common(): print('  %-28s %5.1f%%'%(k,100*v/ex))
# mutual pairs
mut=[]
for (a,b),v in pair.items():
    if a<b and (b,a) in pair and a!=b:
        mut.append((min(v,pair[(b,a)]),a,b,v,pair[(b,a)]))
mut.sort(reverse=True)
mt=sum(v+w for _,a,b,v,w in mut)
print('  mutual cross-copying pairs: %d, %.1f%% of copies; top:'%(len(mut),100*mt/ex))
for m,a,b,v,w in mut[:4]: print('    %s <-> %s  %d / %d   self %d / %d'%(seqs[a],seqs[b],v,w,pair.get((a,a),0),pair.get((b,b),0)))
print('  top copy relations:')
for (a,b),v in pair.most_common(8): print('    %5.1f%%  %s  copies  %s'%(100*v/ex,seqs[a],'ITSELF' if a==b else seqs[b]))
# who copies passengers
pc=collections.Counter()
for (a,b),v in pair.items():
    if not cop(b): pc[a]+=v
ps=sum(pc.values())
print('  main hosts of passengers:')
for a,v in pc.most_common(4): print('    %5.1f%%  %s'%(100*v/ps,seqs[a]))
