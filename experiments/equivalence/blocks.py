import os, sys, collections
def blocks(p):
    b=[]; cur=None
    for l in open(p, errors="replace"):
        if l.startswith("PAIR"):
            cur=[l]; b.append(cur)
        else: cur.append(l)
    return b
kinds=collections.Counter(); ex={}
tot=0; ndiff=0
for d in sorted(os.listdir("work")):
    A=blocks(f"work/{d}/cpp.txt"); B=blocks(f"work/{d}/smol.txt")
    for a,b in zip(A,B):
        tot+=1
        if a==b: continue
        ndiff+=1
        i=next(k for k in range(min(len(a),len(b))) if a[k]!=b[k]) if any(x!=y for x,y in zip(a,b)) else min(len(a),len(b))
        la=a[i] if i<len(a) else "<none>"; lb=b[i] if i<len(b) else "<none>"
        if la.startswith("END") and lb.startswith("END"):
            sa=la.split(); sb=lb.split()
            k="END:"+("steps" if sa[1]!=sb[1] else "act" if sa[-2]!=sb[-2] else "pas")
        elif la.startswith("PAIR"): k="PAIR(bind)"
        else: k="TRACE"
        kinds[k]+=1
        ex.setdefault(k,[]).append((d,a[0].strip(),i,la.strip()[:160],lb.strip()[:160], max(len(x) for x in a[-1].split()) ))
print("pairs in diff files",tot,"differing",ndiff,dict(kinds))
for k,v in ex.items():
    print("==",k,len(v))
    for e in v[:4]: print("  ",e)
