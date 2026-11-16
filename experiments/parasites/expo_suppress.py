import sys, re, glob, gzip, collections, subprocess, tempfile, os
R=os.path.expanduser("~/smol-src/target/release/smol-react")
def probe(seqs, flags):
    with tempfile.NamedTemporaryFile("w",suffix=".txt",delete=False) as f:
        for i,s in enumerate(seqs): f.write(f"R m{i} {s}\n")
        n=f.name
    o=subprocess.run([R,n,"5000"]+flags,capture_output=True,text=True).stdout; os.unlink(n); return o
for d in sorted(glob.glob(os.path.expanduser("~/smol/runs/suppress/*/s*"))):
    for step in (300000,1000000,2000000):
        ck=f"{d}/checkpoint_{step}.txt.gz"
        if not os.path.exists(ck): continue
        c=collections.Counter()
        with gzip.open(ck,"rt") as f:
            for l in f:
                if l.startswith("MOLS"): break
            for l in f:
                p=l.split()
                if p[5]=="F": c[p[4].replace(".","")]+=1
        top=[s for s,_ in c.most_common(12)]
        out=probe(top,["--self"]); S=None
        for i,s in enumerate(top):
            m=re.search(rf"^m{i}\s+initiates with m{i}\s.*\n(?:\s+->.*\n)*",out,re.M)
            if m and re.search(rf"-> product m{i}\b",m.group(0)): S=s; break
        if not S: continue
        tr=probe([S],["--steptrace","--self"])
        nc=sum(1 for l in tr.splitlines() if l.startswith("T ") and l.split()[2]=="=")
        pl=sum(len(l.split()[2]) for l in tr.splitlines() if l.startswith("PROD"))
        print(d.split("/")[-2],d.split("/")[-1],step,len(S),f"{nc/pl:.2f}" if pl else "nan",S[:40],flush=True)
