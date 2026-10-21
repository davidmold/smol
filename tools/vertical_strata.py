import struct,zlib,sys,collections
for d in sys.argv[1:]:
    b=open(d+'/live.grid','rb').read()
    step,pw,ph,zl=struct.unpack('<4I',b[:16]); a=struct.unpack('<%dI'%(pw*ph),zlib.decompress(b[16:16+zl]))
    gy=12; gz=(ph+1)//(gy+1)
    seqs={}
    for l in open(d+'/species.tsv'):
        p=l.rstrip('\n').split('\t')
        if p[0].isdigit(): seqs[int(p[0])]=p[-1]
    w={}
    for l in open(d+'/decay_weights.tsv'):
        p=l.split()
        if p[0]!='symbol': w[p[0]]=float(p[1])
    f=lambda s: sum(w.get(c,1) for c in s)/max(1,len(s))
    dep=collections.defaultdict(list); cop=collections.defaultdict(list); occ=[]; tops=collections.Counter()
    for z in range(gz):
        n=0
        for x in range(pw):
            col=[a[(z*(gy+1)+y)*pw+x]&0x7fffffff for y in range(gy) if a[(z*(gy+1)+y)*pw+x]]
            n+=len(col)
            for k,v in enumerate(col):
                s=seqs.get(v)
                if s is None: continue
                if k==0: tops[s]+=1
                dep[min(k,4)].append(f(s)); cop[min(k,4)].append('=' in s)
        occ.append(n)
    allv=[v for l in dep.values() for v in l]
    print(d,'step',step,'per slice',occ,'mean dfac %.2f'%(sum(allv)/max(1,len(allv))),'depth:',' '.join('%d:%.2f/%.0f%%'%(k,sum(v)/len(v),100*sum(cop[k])/len(v)) for k,v in sorted(dep.items())),'| tops',tops.most_common(3))
