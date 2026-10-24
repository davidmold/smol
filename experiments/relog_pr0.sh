#!/bin/bash
# Rerun the 24 surviving pr0 seeds to 1M with OFFSPRINGLOG (same seeds and chemistry => same runs).
R=~/smol/runs; BIN=~/smol/target/release/smol
while [ ! -e $R/dep/DONE ]; do sleep 60; done
J=$R/pr0log_jobs.txt; : > $J
for f in $R/pr0/seed*/log.txt; do
  grep -q "Step 1000000 done" $f || continue
  s=$(basename $(dirname $f)); d=$R/pr0log/$s; mkdir -p $d
  { grep -vE "^(OFFSPRINGLOG|AGENT|GRIDPOS)" $(dirname $f)/run.conf; echo "OFFSPRINGLOG 1"; grep -E "^(AGENT|GRIDPOS)" $(dirname $f)/run.conf; } > $d/run.conf
  echo $d >> $J
done
cat $J | xargs -P 10 -I{} sh -c 'cd {} && [ -s log.txt ] || '$BIN' run run.conf > log.txt 2>&1'
echo DONE > $R/pr0log/DONE
