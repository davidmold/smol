#!/bin/bash
# prep.sh STEP: convert the 8 published runs at STEP and run a 20k-step classification window
set -e
STEP=$1; B=~/smol-src/target/release/smol; T=~/smol-src/tools; O=~/smol/runs/parasites/s$STEP
for r in 1 2 3 4 5 6 7 8; do
  d=$O/r$r/src; mkdir -p $d
  python3 $T/pub_to_ckpt.py ~/published/bio/biopaperdata/raw/run$r/out1_$STEP.conf $d --nsteps $((STEP+20000)) > $d/convert.txt
  ( cd $d && $B resume checkpoint_$STEP.txt.gz --reseed 1 > log.txt 2>&1 ) &
done
wait
