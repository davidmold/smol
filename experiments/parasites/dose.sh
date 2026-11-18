#!/bin/bash
# dose.sh: measures for the graded PARASITEFREE 2 runs (PFPROB = blocking probability q).
# q = 0 is the control runs (ctl), q = 1 is PARASITEFREE 2 (pf2); same seeds 201-210 throughout.
cd ~/smol/runs/pfree
V=~/smol/runs/doseview; mkdir -p $V
ln -sfn ../pfree/ctl $V/q000; ln -sfn ../pfree/pfq25 $V/q025; ln -sfn ../pfree/pfq50 $V/q050
ln -sfn ../pfree/pfq75 $V/q075; ln -sfn ../pfree/pf2 $V/q100
cd ~/smol/runs
python3 complexity.py doseview > pfree/dose_complexity.txt 2>&1
python3 expo_dose.py "~/smol/runs/doseview/q*/s*" > pfree/dose_expo.txt 2>&1
python3 scanpath.py doseview/q*/s2* > pfree/dose_scanpath.txt 2>&1
for d in doseview/q*/s*; do echo "$d $(tail -1 $d/log.txt | grep -oE 'EXTINCT at step [0-9]+|^Step [0-9]+' )"; done > pfree/dose_survival.txt
for q in q025 q050 q075 q100; do for d in doseview/$q/s*; do
  awk -v d=$d '/^Step/ {for(i=1;i<=NF;i++){if($i=="blocked")b+=$(i+1); if($i=="pfpassed")p+=$(i+1)}} END{print d, b, p}' $d/log.txt
done; done > pfree/dose_load.txt
touch pfree/dose.done
