import sys,statistics,collections,math
for fn in sys.argv[1:]:
    T=collections.defaultdict(list)
    for l in open(fn):
        p=l.split()
        if len(p)==9 and p[0].isdigit():
            v=[float(x) for x in p]; T[int(p[0])].append(v)
    print('==',fn,'(median across runs; parasitic and mutual are shares of productive cross-pairs)')
    print('%8s %4s %6s %7s %6s %9s %6s %6s'%('step','runs','N','selfrep','mutual','parasitic','rtime','toplen'))
    for t in sorted(T):
        r=T[t]; m=lambda k: statistics.median([x[k] for x in r if not math.isnan(x[k])] or [float('nan')])
        print('%8d %4d %6.0f %7.2f %6.2f %9.2f %6.0f %6.0f'%(t,len(r),m(1),m(3),m(4),m(5),m(6),m(7)))
