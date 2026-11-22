#!/bin/bash
# confirm_pipeline.sh: the whole confirmation round (CONFIRM_PLAN.md), unattended on aida.
# 1. runs (confirm_setup.py's jobs, 12 at a time); 2. time-shift matrices for cf_gmt, cf_gmq,
# cf_gmd (200k-2M) and for gmt6, gmq6 (2M-6M every 0.5M); 3. run-level statistics and mechanism.
R=$HOME/smol/runs; E=$HOME/smol-src/experiments/parasites; cd $R
python3 $E/confirm_setup.py > confirm_setup.log 2>&1
xargs -P 12 -I{} bash -c "{}" < confirm_jobs.txt; touch confirm_runs.done
for s in cf_gmt cf_gmq cf_gmd; do
  CLASSIFY=classify_gm.py GMFLAGS=--qsearch bash $E/timeshift_smol.sh $R/$s $R/ts_$s > ts_$s.log 2>&1
  python3 $E/gm_summary.py $s ts_$s > ${s}_summary.txt 2>&1
  python3 $E/timeshift_agg.py ts_$s 10000 > ts_$s.results.txt 2>&1
  python3 $E/timeshift_extra.py ts_$s > ts_$s.extra.txt 2>&1
  python3 $E/gmt_mech.py ts_$s > ${s}_mech.txt 2>&1
done
python3 $E/confirm_stats.py ts_cf_gmt ts_cf_gmq ts_cf_gmd > confirm_primary.txt 2>&1
for s in gmt6 gmq6; do
  TIMES="2000000 2500000 3000000 3500000 4000000 4500000 5000000 5500000 6000000" \
    CLASSIFY=classify_gm.py GMFLAGS=--qsearch bash $E/timeshift_smol.sh $R/$s $R/ts_$s > ts_$s.log 2>&1
  python3 $E/timeshift_agg.py ts_$s 10000 > ts_$s.results.txt 2>&1
done
python3 $E/confirm_stats.py ts_gmt6 ts_gmq6 > confirm_long.txt 2>&1
df -h ~ | tail -1 > confirm_pipeline.done
