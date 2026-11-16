#!/bin/bash
# pswap.sh NSTEPS: parasite swaps (the other half of the 2x2): 1M worlds get own vs 2M parasites,
# 2M worlds get own vs 1M parasites; 5 seeds, all 8 runs. Waits for the replicator swaps.
N=$1; B=~/smol-src/target/release/smol; P=~/smol/runs/parasites; O=$P/swap
while [ ! -f $O/swap.done ]; do sleep 60; done
jobs=$O/pjobs.txt; : > $jobs
for r in 1 2 3 4 5 6 7 8; do
  a=$P/s1000000/r$r/src; b=$P/s2000000/r$r/src
  for seed in 1 2 3 4 5; do
    for spec in "pself1:$a:1000000:$a:1000000" "plate2:$a:1000000:$b:2000000" "pself2:$b:2000000:$b:2000000" "pearly1:$b:2000000:$a:1000000"; do
      IFS=: read t rs rstep ds dstep <<< "$spec"
      d=$O/r$r/$t/s$seed
      python3 ~/smol-src/experiments/parasites/swap.py $rs $rstep $ds $dstep $d $seed $N P
      echo "cd $d && $B resume checkpoint_$rstep.txt.gz --reseed $((seed+10)) > log.txt 2>&1" >> $jobs
    done
  done
done
xargs -P 11 -I{} bash -c "{}" < $jobs
touch $O/pswap.done
