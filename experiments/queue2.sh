#!/bin/bash
# Queue after partnership batch 1 (waits for runs/partners/DONE). Runs 9 at a time on aida.
#  Q1 replay/   strand0 seeds where borrower partnerships established: rerun from scratch with
#               a permanent checkpoint every 10k steps and OFFSPRINGLOG, to 100k past
#               establishment. Same seed + same chemistry => same run as the original.
#  Q2 pub/      matched to Hickinbotham, Stepney & Hogeweg (2021): 12,500 cells as 125x100 or
#               250x50, 100 seed replicases as a 10x10 block or 50 adjacent pairs at random
#               places, old partner and turn rules, 2M steps, 5 seeds per treatment.
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

# Q2
SEEDREP='WWGEWLHHHRLUEUWJJJRJXUUUDYGRHJLRWWRE$BLUBO^B>C$=?>$$BLUBO%}OYHOB'
for grid in 125x100 250x50; do for place in block pairs; do for s in 1101 2202 3303 4404 5505; do
  d=$R/pub/${grid}_${place}/seed$s; mkdir -p $d
  python3 - "$d/run.conf" "$grid" "$place" "$s" "$SEEDREP" <<'EOF'
import sys, random
out, grid, place, seed, rep = sys.argv[1:]
gx, gy = map(int, grid.split("x")); rng = random.Random(int(seed))
L = [f"GRIDX {gx}", f"GRIDY {gy}", "DECAY 0.0005", "MUTATE 0.0002", f"RANDSEED {seed}",
     "NSTEPS 2000000", "COMPAT_PARTNER 1", "ENCOUNTER_TURN 1", "OFFSPRINGLOG 1",
     "POPDYEVERY 1000", "SNAPEVERY 5000", "COPYEVERY 1000", "CKPTEVERY 100000",
     "KEEPCKPT 2", "CKPTKEEPEVERY 100000"]
cells = []
if place == "block":
    cx, cy = gx // 2, gy // 2
    cells = [(cx - 5 + i, cy - 5 + j) for i in range(10) for j in range(10)]
else:
    taken = set()
    while len(cells) < 100:
        x, y = rng.randrange(gx), rng.randrange(gy); x2 = (x + 1) % gx
        if (x, y) in taken or (x2, y) in taken:
            continue
        taken |= {(x, y), (x2, y)}; cells += [(x, y), (x2, y)]
L.append(f"AGENT  {rep} 100 Q"); L += [f"GRIDPOS {x} {y}" for x, y in cells]
open(out, "w").write("\n".join(L) + "\n")
EOF
  echo "run $d" >> $J
done; done; done

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
