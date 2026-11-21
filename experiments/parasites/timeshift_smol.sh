#!/bin/bash
# timeshift_smol.sh SETDIR OUT: time-shift matrix for smol's own runs (SETDIR/sN with checkpoints
# every 100k). For each run and each time t in 200k..2M where the run is alive: resume a 20k
# classification window from checkpoint t (classify.py), then for every host time H and parasite
# time P replace H's parasites with P's (swap.py ... P), 2 seeds, 20k steps. Same measure as
# timeshift.sh; summarise with timeshift_agg.py OUT.
set -e
SET=$1; O=$2; mkdir -p $O
B=~/smol-src/target/release/smol; T=~/smol-src/tools; E=~/smol-src/experiments/parasites
TIMES="200000 400000 600000 800000 1000000 1200000 1400000 1600000 1800000 2000000"
N=20000
# 1. classification windows
wjobs=$O/wjobs.txt; : > $wjobs
for run in $SET/s*; do
  r=$(basename $run)
  for t in $TIMES; do
    [ -f $run/checkpoint_$t.txt.gz ] || continue
    src=$O/$r/w$((t/1000))k
    [ -f $src/classes.tsv ] && continue
    echo "python3 $T/ckpt_edit.py $run checkpoint_$t.txt.gz $src --set NSTEPS $((t+20000)) --set POPDYEVERY 1000 --set COPYEVERY 1000 > /dev/null 2>&1 && cd $src && $B resume checkpoint_$t.txt.gz --reseed 1 > log.txt 2>&1; python3 $E/${CLASSIFY:-classify.py} $src $t > /dev/null" >> $wjobs
  done
done
xargs -P 12 -I{} bash -c "{}" < $wjobs || true
touch $O/windows.done
# Each swap run deletes its copies of species.tsv, copies.tsv and grid.snap when done: with
# smol's own 2M runs the species list is large and 1,000+ copies filled the disk (21 Nov).
# 2. swaps (only worlds that have a classes.tsv and some parasites)
jobs=$O/jobs.txt; : > $jobs
for rd in $O/s*; do
  r=$(basename $rd)
  for h in $TIMES; do
    hs=$rd/w$((h/1000))k; [ -s $hs/classes.tsv ] || continue
    grep -q " P " <(awk '{print " "$2" "}' $hs/classes.tsv) || continue
    for p in $TIMES; do
      ps=$rd/w$((p/1000))k; [ -s $ps/classes.tsv ] || continue
      grep -q " P " <(awk '{print " "$2" "}' $ps/classes.tsv) || continue
      for seed in 1 2; do
        d=$O/$r/h$((h/1000))k_p$((p/1000))k/s$seed
        echo "python3 $E/swap.py $hs $h $ps $p $d $seed $N P > /dev/null && cd $d && $B resume checkpoint_$h.txt.gz --reseed $((seed+10)) > log.txt 2>&1; rm -f species.tsv copies.tsv grid.snap" >> $jobs
      done
    done
  done
done
wc -l $jobs
xargs -P 12 -I{} bash -c "{}" < $jobs || true
touch $O/all.done
