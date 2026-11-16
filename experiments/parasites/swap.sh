#!/bin/bash
# swap.sh NSTEPS: forward (1M recipients get own vs 2M replicators) and reverse (2M recipients
# get own vs 1M replicators), 5 seeds each, all 8 runs
N=$1; B=~/smol-src/target/release/smol; P=~/smol/runs/parasites; O=$P/swap; mkdir -p $O; jobs=$O/jobs.txt; : > $jobs
for r in 1 2 3 4 5 6 7 8; do
  a=$P/s1000000/r$r/src; b=$P/s2000000/r$r/src
  python3 ~/smol-src/experiments/parasites/classify.py $b 2000000 >> $O/classes_2M.txt
  for seed in 1 2 3 4 5; do
    for spec in "self1:$a:1000000:$a:1000000" "late2:$a:1000000:$b:2000000" "self2:$b:2000000:$b:2000000" "early1:$b:2000000:$a:1000000"; do
      IFS=: read t rs rstep ds dstep <<< "$spec"
      d=$O/r$r/$t/s$seed
      python3 ~/smol-src/experiments/parasites/swap.py $rs $rstep $ds $dstep $d $seed $N
      echo "cd $d && $B resume checkpoint_$rstep.txt.gz --reseed $((seed+10)) > log.txt 2>&1" >> $jobs
    done
  done
done
xargs -P 11 -I{} bash -c "{}" < $jobs
touch $O/swap.done
