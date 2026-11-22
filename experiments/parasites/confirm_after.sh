#!/bin/bash
# confirm_after.sh: after confirm_pipeline.sh, the widened mechanism assay (addendum point 2) on
# the new sets and the extensions, then the run-level statistics again with the completeness check.
R=$HOME/smol/runs; E=$HOME/smol-src/experiments/parasites; cd $R
while [ ! -f confirm_pipeline.done ]; do sleep 300; done
~/smenv/bin/python $E/mech2.py ts_cf_gmt ts_cf_gmq ts_cf_gmd -j 12 > mech2_confirm.txt 2>&1
~/smenv/bin/python $E/mech2.py ts_gmt6 ts_gmq6 -j 12 > mech2_long.txt 2>&1
python3 $E/confirm_stats.py ts_cf_gmt ts_cf_gmq ts_cf_gmd > confirm_primary.txt 2>&1
python3 $E/confirm_stats.py ts_cf_gmd > confirm_gmd_alone.txt 2>&1
python3 $E/confirm_stats.py ts_gmt6 ts_gmq6 > confirm_long.txt 2>&1
touch confirm_after.done
