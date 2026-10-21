import sys,collections,struct,zlib,json
def stats(f,every=100000):
    out={};cur=None;d=collections.Counter()
    for line in open(f):
        p=line.strip().split(',')
        if len(p)<3: continue
        s=int(p[0])
        if s!=cur:
            if cur is not None and cur%every==0: out[cur]=d
            cur=s;d=collections.Counter()
        d[int(p[1])]+=int(p[2])
    return out
def summ(d):
    n=sum(d.values()); top=d.most_common(1)[0][1]
    simp=sum((c/n)**2 for c in d.values())
    return n,len(d),round(top/n,3),round(1/simp,1)
for run in ('strand_400x300','gradient_400x300','classic_400x300'):
    D='/home/davidmold/smol/runs/big/%s/'%run
    seqs={}
    for line in open(D+'species.tsv'):
        p=line.rstrip('\n').split('\t')
        if p[0].isdigit(): seqs[int(p[0])]=p[-1]
    print('==',run)
    st=stats(D+'popdy001.dat')
    for s in sorted(st):
        if s%200000==0 or s==max(st): 
            d=st[s]; n=sum(d.values())
            ml=sum(len(seqs.get(k,''))*v for k,v in d.items())/n
            print(s,*summ(d),'meanlen',round(ml,1))
    b=open(D+'live.grid','rb').read()
    step,gx,gy,zl=struct.unpack('<4I',b[:16])
    a=struct.unpack('<%dI'%(gx*gy),zlib.decompress(b[16:16+zl]))
    print('grid step',step)
    # bands by distance from x=200 (warmth)
    for lo,hi in ((0,20),(20,60),(60,100),(100,140),(140,180),(180,200)):
        c=collections.Counter();cells=0
        for y in range(gy):
            for x in range(gx):
                dist=abs(x-200)
                if lo<=dist<hi:
                    cells+=1;v=a[y*gx+x]
                    if v: c[v&0x7fffffff]+=1
        n=sum(c.values()); ml=sum(len(seqs.get(k,''))*v for k,v in c.items())/max(n,1)
        print(' dist',lo,hi,'occ',round(n/cells,2),'sp',len(c),'meanlen',round(ml,1),[(round(v/n,2),seqs.get(k,'?')[:22]) for k,v in c.most_common(3)])
