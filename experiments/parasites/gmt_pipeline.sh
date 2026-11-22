#!/bin/bash
# gmt_pipeline.sh: after the tail-check (gmt) and QSEARCH-control (gmq) gene-and-machine runs,
# time-shift matrices for both (gene-parasites swapped; classify_gm with --qsearch), summaries.
R=$HOME/smol/runs; E=$HOME/smol-src/experiments/parasites
while [ ! -f $R/gmt/all.done ] || [ ! -f $R/gmt_runs.done ]; do sleep 60; done
cd $R
for s in gmt gmq; do
  CLASSIFY=classify_gm.py GMFLAGS=--qsearch bash $E/timeshift_smol.sh $R/$s $R/ts_$s > $R/ts_$s.log 2>&1
  python3 $E/gm_summary.py $s ts_$s > ${s}_summary.txt 2>&1
  python3 timeshift_agg.py ts_$s 10000 > ts_$s.results.txt 2>&1
  python3 timeshift_extra.py ts_$s > ts_$s.extra.txt 2>&1
done
df -h ~ | tail -1 > gmt_pipeline.done
