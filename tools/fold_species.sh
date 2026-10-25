#!/bin/bash
# fold a plain species.tsv into species.tsv.gz (as the first gzip member), for finished runs
d=$1; cd $d || exit 0
[ -f species.tsv ] || exit 0
if [ -f species.tsv.gz ]; then
  nice gzip -c species.tsv > species.tsv.new.gz && cat species.tsv.gz >> species.tsv.new.gz && \
  mv species.tsv.new.gz species.tsv.gz && rm species.tsv
else
  nice gzip species.tsv
fi
