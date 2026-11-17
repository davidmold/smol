import re, subprocess, tempfile, os
R=os.path.expanduser("~/smol-src/target/release/smol-react")
for l in open(os.path.expanduser("~/smol/runs/compete/pairs.tsv")):
    run,step,A,D,par=l.rstrip("\n").split("\t")[:5]
    ps=[p for p in par.split(",") if p]
    with tempfile.NamedTemporaryFile("w",suffix=".txt",delete=False) as f:
        f.write(f"R A {A}\nR D {D}\n"); [f.write(f"R p{i} {p}\n") for i,p in enumerate(ps)]; n=f.name
    out=subprocess.run([R,n,"5000"],capture_output=True,text=True).stdout; os.unlink(n)
    res={}; cur=None
    for x in out.splitlines():
        m=re.match(r"(\S+)\s+initiates with (\S+)\s+P\(bind\)=([\d.]+).*steps=(\d+)",x)
        if m: cur=(m.group(1),m.group(2)); res[cur]=[float(m.group(3)),set(),int(m.group(4))]; continue
        m=re.match(r"\s+-> (?:product|passive made) (\S+)",x)
        if m and cur: res[cur][1].add(m.group(1))
    def score(rep):
        s=0
        for i in range(len(ps)):
            pi=f"p{i}"; e=max((res.get((rep,pi),[0,set()]),res.get((pi,rep),[0,set()])),key=lambda v:v[0])
            cop = pi in res.get((rep,pi),[0,set()])[1] or pi in res.get((pi,rep),[0,set()])[1]
            s+= e[0]*cop
        return s/len(ps) if ps else 0
    def tsteps(rep):
        v=[]
        for i in range(len(ps)):
            pi=f"p{i}"
            for k in ((rep,pi),(pi,rep)):
                if k in res and res[k][0]>0: v.append(res[k][2])
        return sum(v)/len(v) if v else 0
    helpA=sum(1 for i in range(len(ps)) for k in (("A",f"p{i}"),(f"p{i}","A")) if "A" in res.get(k,[0,set()])[1])
    helpD=sum(1 for i in range(len(ps)) for k in (("D",f"p{i}"),(f"p{i}","D")) if "D" in res.get(k,[0,set()])[1])
    print(run, "reactions where the replicator itself is copied: A %d D %d;"%(helpA,helpD), end=" ")
    print(run, "mean steps bound to parasites: A %.0f D %.0f;"%(tsteps("A"),tsteps("D")), "parasite exposure (mean P(bind) x copied): A %.2f  D %.2f"%(score("A"),score("D")), " AxD:", sorted(res.get(("A","D"),[0,set()])[1]|res.get(("D","A"),[0,set()])[1]))
