#!/usr/bin/env python3
"""confirm_pshare.py: exploratory (not in the plan): gene-parasite molecule share at 0.2M and 2M per run in the confirmation sets, change within sets and 2M share between sets."""
import glob,os,sys
from scipy import stats
R=os.path.expanduser("~/smol/runs")
def share(ts,t):
    out=[]
    for c in sorted(glob.glob(f"{R}/{ts}/s*/w{t}k/classes.tsv")):
        n={}
        for l in open(c):
            a,b,k=l.split(); n[b]=n.get(b,0)+int(k)
        tot=sum(n.values()); out.append(n.get("P",0)/tot)
    return out
d={}
for s in ("cf_gmt","cf_gmq","cf_gmd"):
    a=share("ts_"+s,200); b=share("ts_"+s,2000); d[s]=(a,b)
    ch=[y-x for x,y in zip(a,b)]
    import statistics as S
    print(s, "P share 0.2M %.2f 2M %.2f  change median %+.3f  falls in %d/%d  Wilcoxon p=%.2g"%(S.median(a),S.median(b),S.median(ch),sum(c<0 for c in ch),len(ch),stats.wilcoxon(ch).pvalue))
for o in ("cf_gmq","cf_gmd"):
    print("2M P share cf_gmt vs",o,"MW p=%.2g"%stats.mannwhitneyu(d["cf_gmt"][1],d[o][1]).pvalue)
