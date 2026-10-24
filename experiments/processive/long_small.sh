#!/bin/bash
# Small classic worlds through the window where the big world's transition happened.
# Waits for the transplant batch, then runs 9 at a time:
#   pr0xx/seedS  the 18 baseline survivors at 2.5M continued to 10M (survivor-conditioned)
#   long/seedS   20 fresh baseline runs from step 0 to 10M (extinction is a competing outcome)
# All with OFFSPRINGLOG and GZIPLOGS, permanent checkpoints every 500k.
R=~/smol/runs; BIN=~/smol/target/release/smol; ED="python3 $HOME/smol/tools/ckpt_edit.py"
while [ ! -e $R/transplant/DONE ]; do sleep 120; done
J=$R/long_jobs.txt; : > $J
for d in $R/pr0x/seed*; do
  [ -e $d/checkpoint_2500000.txt.gz ] || continue
  o=$R/pr0xx/$(basename $d)
  [ -e $o/checkpoint_2500000.txt.gz ] || $ED $d checkpoint_2500000.txt.gz $o --set NSTEPS 10000000 --set GZIPLOGS 1 \
     --set CKPTEVERY 100000 --set KEEPCKPT 3 --set CKPTKEEPEVERY 500000 > /dev/null
  echo "resume $o checkpoint_2500000.txt.gz" >> $J
done
for s in 58058 59059 60060 61061 62062 63063 64064 65065 66066 67067 68068 69069 70070 71071 72072 73073 74074 75075 76076 77077; do
  o=$R/long/seed$s; mkdir -p $o
  { grep -vE "^(RANDSEED|NSTEPS|AGENT|GRIDPOS)" $R/pr0/seed18018/run.conf
    printf "RANDSEED $s\nNSTEPS 10000000\nOFFSPRINGLOG 1\nGZIPLOGS 1\nCKPTEVERY 100000\nKEEPCKPT 3\nCKPTKEEPEVERY 500000\n"
    grep -E "^(AGENT|GRIDPOS)" $R/pr0/seed18018/run.conf; } > $o/run.conf
  echo "run $o run.conf" >> $J
done
go() { set -- $1; cd $2; [ -s log.txt ] && return 0
  if [ $1 = run ]; then $BIN run $3 > log.txt 2>&1; else $BIN resume $3 > log.txt 2>&1; fi; }
export -f go; export BIN
cat $J | xargs -P 9 -I{} bash -c 'go "{}"'
echo DONE > $R/long_DONE
