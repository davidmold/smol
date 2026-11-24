#!/bin/bash
# pw_pipeline.sh: find password transitions, set up and run the contests, score (pw_compete.py).
R=$HOME/smol/runs; E=$HOME/smol-src/experiments/parasites; mkdir -p $R/pwc; cd $R/pwc
python3 $E/pw_compete.py find > find.txt 2>&1 || exit 1
python3 $E/pw_compete.py setup > setup.txt 2>&1 || exit 1
xargs -P 12 -I{} bash -c "{}" < jobs.txt
~/smenv/bin/python $E/pw_compete.py score > score.txt 2>&1
touch pwc.done
