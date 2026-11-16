#!/bin/bash
# treat.sh STEP NSTEPS: build the edited worlds and run them (base, noP, noN, cullP, cullN x 3 seeds)
STEP=$1; N=$2; B=~/smol-src/target/release/smol; T=~/smol-src/tools; O=~/smol/runs/parasites/s$STEP
jobs=$O/jobs.txt; : > $jobs
for r in 1 2 3 4 5 6 7 8; do
  src=$O/r$r/src; ck=$src/checkpoint_$STEP.txt.gz
  python3 ~/smol-src/experiments/parasites/classify.py $src $STEP >> $O/classes_summary.txt
  tot=$(awk '{t+=$3} END{print t}' $src/classes.tsv)
  for c in P N; do
    ids=$(awk -v c=$c '$2==c{printf " --remove-species %s",$1}' $src/classes.tsv)
    n=$(awk -v c=$c '$2==c{t+=$3} END{print t+0}' $src/classes.tsv)
    keep=$(python3 -c "print(1-$n/$tot)")
    for seed in 2 3 4; do
      for t in no$c cull$c; do
        d=$O/r$r/$t/s$seed; mkdir -p $d
        if [ $t = no$c ]; then ops="$ids"; else ops="--cull $keep"; fi
        python3 $T/ckpt_edit.py $src checkpoint_$STEP.txt.gz $d $ops --set NSTEPS $((STEP+N)) --set POPDYEVERY 1000 --rng-seed $seed > $d/edit.txt 2>&1
        echo "cd $d && $B resume checkpoint_$STEP.txt.gz --reseed $seed > log.txt 2>&1" >> $jobs
      done
    done
  done
  for seed in 2 3 4; do
    d=$O/r$r/base/s$seed; mkdir -p $d
    python3 $T/ckpt_edit.py $src checkpoint_$STEP.txt.gz $d --set NSTEPS $((STEP+N)) --set POPDYEVERY 1000 --rng-seed $seed > $d/edit.txt 2>&1
    echo "cd $d && $B resume checkpoint_$STEP.txt.gz --reseed $seed > log.txt 2>&1" >> $jobs
  done
done
sort -u $jobs -o $jobs
xargs -P 10 -I{} bash -c "{}" < $jobs
touch $O/treat.done
