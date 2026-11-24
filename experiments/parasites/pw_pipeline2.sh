#!/bin/bash
# pw_pipeline2.sh: second round of password contests (pw_compete2.py); rebuild and setup done first.
cd $HOME/smol/runs/pwc2 && xargs -P 12 -I{} bash -c "{}" < jobs.txt
~/smenv/bin/python $HOME/smol-src/experiments/parasites/pw_compete2.py score > score.txt 2>&1
touch pwc2.done
