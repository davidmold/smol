#!/bin/bash
# pw_replay.sh: lineage replays (pw_replay.py): setup, 12 at a time, score.
cd $HOME/smol/runs && mkdir -p pwr && cd pwr
python3 $HOME/smol-src/experiments/parasites/pw_replay.py setup > setup.txt 2>&1 || exit 1
xargs -P 12 -I{} bash -c "{}" < jobs.txt
~/smenv/bin/python $HOME/smol-src/experiments/parasites/pw_replay.py score > score.txt 2>&1
touch pwr.done
