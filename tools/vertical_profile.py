import struct,zlib,sys,collections
d=sys.argv[1]
b=open(d+'/live.grid','rb').read()
step,gx,gy,zl=struct.unpack('<4I',b[:16])
a=struct.unpack('<%dI'%(gx*gy),zlib.decompress(b[16:16+zl]))
seqs={}
for l in open(d+'/species.tsv'):
    p=l.rstrip('\n').split('\t')
    if p[0].isdigit(): seqs[int(p[0])]=p[-1]
w={}
try:
    for l in open(d+'/decay_weights.tsv'):
        p=l.split()
        if len(p)==2 and p[0]!='symbol': w[p[0]]=float(p[1])
except: pass
print('step',step)
for y in range(gy):
    row=[a[y*gx+x]&0x7fffffff for x in range(gx)]
    n=sum(1 for v in row if v)
    if n: 
        c=collections.Counter(v for v in row if v)
        s=seqs.get(c.most_common(1)[0][0],'')
        f=sum(w.get(ch,1) for ch in s)/max(1,len(s))
        print('h%2d occ %3d top %-22s dfac %.2f'%(gy-1-y,n,s[:22],f))
