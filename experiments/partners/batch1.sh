#!/bin/bash
# Partnership isolation, batch 1: each pair alone, pair vs the plain copier, pair with a
# parasite; mutation off and on; 10 seeds; 200k steps on 125x100. Members alone: 2 seeds.
# Usage on aida: bash batch1.sh  (writes ~/smol/runs/partners/, runs 9 at a time)
set -e
R=~/smol/runs/partners; M=~/smol/experiments/partners/mkconf.py; BIN=~/smol/target/release/smol
COP='RERERERE^BQQ$=?>%}'
declare -A A B P
A[s1]='S^BO$=J?>G}';      B[s1]='F^O>';      P[s1]='FIO>'
A[s2]='RER$Y$^BQQ$=?>H}'; B[s2]='ERE^AEBP>'; P[s2]='ERE^'
A[s3]='FSF^BQ$=?>G}';     B[s3]='SES^AQ';    P[s3]='SES'
A[g1]='FUS^YF^BPQ$=?L>}'; B[g1]='SIF^ABQ';   P[g1]='SIF'
mkdir -p $R; : > $R/jobs.txt
add() { local d=$R/$1/seed$2; mkdir -p $d; shift 2; python3 $M $d/run.conf "$@"; echo $d >> $R/jobs.txt; }
for p in s1 s2 s3 g1; do
  for s in 18018 19019; do
    add ${p}_Aalone $s --seed $s "${A[$p]}@62,50"
    add ${p}_Balone $s --seed $s "${B[$p]}@62,50"
  done
  for m in 0 0.0002; do
    for s in 18018 19019 20020 21021 22022 23023 24024 25025 26026 27027; do
      add ${p}_AB_m$m $s --seed $s --mutate $m "${A[$p]}+${B[$p]}@62,50"
      add ${p}_ABcop_m$m $s --seed $s --mutate $m "${A[$p]}+${B[$p]}@35,50" "$COP@90,50"
      add ${p}_ABpar_m$m $s --seed $s --mutate $m "${A[$p]}+${B[$p]}+${P[$p]}@62,50"
    done
  done
done
cat $R/jobs.txt | xargs -P 9 -I{} sh -c 'cd {} && [ -e log.txt ] || '$BIN' run run.conf > log.txt 2>&1'
echo DONE > $R/DONE
