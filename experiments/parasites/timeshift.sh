#!/bin/bash
# timeshift.sh: full time-shift matrix of parasites against host worlds, all 8 published runs.
# For each run, host time H and parasite time P in 200k..2M (step 200k): the H world's parasites
# (class P, classify.py over a 20k window) are replaced by parasites drawn by abundance from the
# P world of the same run (swap.py ... P). H = P is the control (same disturbance, own parasites).
# 3 seeds, 20k steps. Lower host population at +10k = the parasites did better.
set -e
B=~/smol-src/target/release/smol; P=~/smol/runs/parasites; O=$P/timeshift; mkdir -p $O
TIMES="200000 400000 600000 800000 1000000 1200000 1400000 1600000 1800000 2000000"
N=20000
# 1. convert and classify every time not done yet (prep.sh runs the 8 runs in parallel)
for t in $TIMES; do
  if [ ! -f $P/s$t/r8/src/classes.tsv ]; then
    bash ~/smol-src/experiments/parasites/prep.sh $t
    for r in 1 2 3 4 5 6 7 8; do
      python3 ~/smol-src/experiments/parasites/classify.py $P/s$t/r$r/src $t > /dev/null
    done
  fi
done
touch $O/prep.done
# 2. one job per cell and seed: build the edited world, then run it
jobs=$O/jobs.txt; : > $jobs
for r in 1 2 3 4 5 6 7 8; do
  for h in $TIMES; do for p in $TIMES; do for seed in 1 2 3; do
    d=$O/r$r/h$((h/1000))k_p$((p/1000))k/s$seed
    echo "python3 ~/smol-src/experiments/parasites/swap.py $P/s$h/r$r/src $h $P/s$p/r$r/src $p $d $seed $N P && cd $d && $B resume checkpoint_$h.txt.gz --reseed $((seed+10)) > log.txt 2>&1" >> $jobs
  done; done; done
done
xargs -P 12 -I{} bash -c "{}" < $jobs || true
touch $O/all.done
