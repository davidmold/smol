#!/bin/bash
# cmpg.sh TAG CPPBIN CONF SMOLFLAGS PAIRFILE
tag=$1; bin=$2; conf=$3; flags=$4; f=$5; n=$(basename $f .txt); w=/tmp/sm_eq/w_$tag/$n; mkdir -p $w /tmp/sm_eq/r_$tag; cd $w
$bin $conf $f 20000 2>/dev/null | grep -aE "^(PAIR|T |PROD|END)" > cpp.txt
~/projects/smol/target/release/smol-react $f 20000 --compat --steptrace $flags > smol.txt
if cmp -s cpp.txt smol.txt; then echo "$n SAME $(grep -c ^PAIR cpp.txt) $(grep -c '^T ' cpp.txt)" > /tmp/sm_eq/r_$tag/$n.txt; cd /tmp/sm_eq; rm -rf $w
else echo "$n DIFF $(grep -c ^PAIR cpp.txt) $(grep -c ^PAIR smol.txt)" > /tmp/sm_eq/r_$tag/$n.txt; fi
