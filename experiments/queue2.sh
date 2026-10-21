#!/bin/bash
# Queue after partnership batch 1 (waits for runs/partners/DONE). Runs 9 at a time on aida.
#  Q1 replay/   strand0 seeds where borrower partnerships established: rerun from scratch with
#               a permanent checkpoint every 10k steps and OFFSPRINGLOG, to 100k past
#               establishment. Same seed + same chemistry => same run as the original.
#  Q2 pub/      the published configurations (Hickinbotham, Stepney & Hogeweg 2021) with the
#               old partner and turn rules, 2M steps, 5 seeds per treatment.
#  Q3 pr0x/     the 24 surviving classic baseline runs (pr0) continued from 1M to 2.5M.
set -e
S=~/smol; R=$S/runs; BIN=$S/target/release/smol; ED="python3 $S/tools/ckpt_edit.py"
while [ ! -e $R/partners/DONE ]; do sleep 60; done
J=$R/queue2_jobs.txt; : > $J

# Q1
declare -A EST=([18018]=110000 [20020]=70000 [23023]=100000 [25025]=140000 [26026]=260000 [24024]=330000 [21021]=450000)
for s in "${!EST[@]}"; do
  d=$R/replay/strand0_seed$s; mkdir -p $d
  { grep -vE "^(NSTEPS|CKPTEVERY|KEEPCKPT|CKPTKEEPEVERY|OFFSPRINGLOG|POPDYEVERY)" $R/strand0/seed$s/run.conf
    echo "NSTEPS $(( ${EST[$s]} + 100000 ))"; echo "CKPTEVERY 10000"; echo "KEEPCKPT 2"
    echo "CKPTKEEPEVERY 10000"; echo "OFFSPRINGLOG 1"
    grep -q "^POPDYEVERY" $R/strand0/seed$s/run.conf && grep "^POPDYEVERY" $R/strand0/seed$s/run.conf || true
  } > $d/run.conf.tmp
  # GRIDPOS lines must stay after their AGENT line: rebuild with the agent block last
  grep -vE "^(AGENT|GRIDPOS)" $d/run.conf.tmp > $d/run.conf
  grep -E "^(AGENT|GRIDPOS)" $R/strand0/seed$s/run.conf >> $d/run.conf; rm $d/run.conf.tmp
  echo "run $d" >> $J
done

# Q2: the published configurations themselves (biopaperdata/script/runN/initrunN.conf):
# box = 125x100 with 121 seeds in a block (runs 1-4), rand = 100x125 with 51 seeds in random
# pairs (runs 5-6), slot = 250x50 with 121 in a block (runs 7-8). Old partner and turn rules,
# 2M steps, 5 random seeds per treatment.
PUB=~/published/bio/biopaperdata/script
for tr in box:1 rand:5 slot:7; do name=${tr%%:*}; src=$PUB/run${tr##*:}/initrun${tr##*:}.conf
  for s in 1101 2202 3303 4404 5505; do
    d=$R/pub/$name/seed$s; mkdir -p $d
    { echo "RANDSEED $s"; echo "NSTEPS 2000000"; echo "COMPAT_PARTNER 1"; echo "ENCOUNTER_TURN 1"
      echo "OFFSPRINGLOG 1"; echo "POPDYEVERY 1000"; echo "SNAPEVERY 5000"; echo "CKPTEVERY 100000"
      echo "KEEPCKPT 2"; echo "CKPTKEEPEVERY 100000"
      grep -E "^(GRIDX|GRIDY|DECAY|MUTATE)" $src
      grep -E "^(AGENT|GRIDPOS)" $src; } > $d/run.conf
    echo "run $d" >> $J
  done
done

# Q3
for f in $R/pr0/seed*/log.txt; do
  grep -q "Step 1000000 done" $f || continue
  src=$(dirname $f); s=$(basename $src); d=$R/pr0x/$s
  [ -e $d/checkpoint_1000000.txt.gz ] || $ED $src checkpoint_1000000.txt.gz $d --set NSTEPS 2500000 --set OFFSPRINGLOG 1 --set CKPTKEEPEVERY 100000 --set CKPTEVERY 100000 > /dev/null
  echo "resume $d" >> $J
done

go() { set -- $1 $2; cd $2; [ -e log.txt ] && return 0
  if [ $1 = run ]; then $BIN run run.conf > log.txt 2>&1; else $BIN resume checkpoint_1000000.txt.gz > log.txt 2>&1; fi; }
export -f go; export BIN
cat $J | xargs -P 9 -I{} bash -c 'go {}'
echo DONE > $R/queue2_DONE
