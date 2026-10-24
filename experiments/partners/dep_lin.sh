#!/bin/bash
cd ~/smol/runs/dep; mkdir -p lin
declare -A A B
A[s1]='S^BO$=J?>G}';      B[s1]='F^O>'
A[s2]='RER$Y$^BQQ$=?>H}'; B[s2]='ERE^AEBP>'
A[s3]='FSF^BQ$=?>G}';     B[s3]='SES^AQ'
A[g1]='FUS^YF^BPQ$=?L>}'; B[g1]='SIF^ABQ'
for p in s1 s2 s3 g1; do for t in ctrl noB noA rescueB; do for d in $p/$t/seed*; do
  [ -s $d/log.txt ] || continue
  python3 ~/smol/experiments/partners/lineage.py $d "A=${A[$p]}" "B=${B[$p]}" --json lin/${p}_${t}_$(basename $d).json > /dev/null 2>&1 || echo "fail $d"
done; done; done
python3 ~/smol/experiments/partners/lineage_agg.py lin > agg.txt
