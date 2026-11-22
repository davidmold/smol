#!/bin/bash
# Archive finished smol run directories from aida to S3 (Glacier Instant Retrieval), via this machine
# so no AWS keys live on aida.
#
#   archive_s3.sh [--delete] RUNPATH...      RUNPATH relative to aida:~/smol/runs, e.g. long/seed59059
#
# Each directory becomes one uncompressed tar (its logs are already gzipped):
#   s3://$SMOL_S3_BUCKET/runs/<RUNPATH>.tar
# The tar stream's SHA-256 and byte count go to s3://$SMOL_S3_BUCKET/manifest/<RUNPATH>.sha256
# and to the local manifest archive_manifest.tsv. The object's size in S3 is checked against the
# byte count; with --delete the directory is removed on aida only after that check passes.
# Restore:  aws s3 cp s3://$SMOL_S3_BUCKET/runs/long/seed59059.tar - | tar -xf -
set -o pipefail
B=${SMOL_S3_BUCKET:?set SMOL_S3_BUCKET to your bucket}
P="--region ${AWS_REGION:-us-east-1}"   # credentials: AWS_PROFILE or the default chain
DEL=0; [ "$1" = "--delete" ] && { DEL=1; shift; }
MAN=${MANIFEST:-$HOME/projects/smol/experiments/archive_manifest.tsv}
tmp=$(mktemp -d)
for r in "$@"; do
  r=${r%/}; key="runs/$r.tar"
  if aws $P s3api head-object --bucket $B --key "$key" >/dev/null 2>&1; then echo "skip $r (already in S3)"; continue; fi
  # in use = a smol process (or a scan) has its working directory or an argument inside the run
  if ssh aida "d=\$(readlink -f ~/smol/runs/$r); for p in \$(pgrep -x smol) \$(pgrep -f transition_scan.py); do c=\$(readlink /proc/\$p/cwd); a=\$(tr '\\0' ' ' < /proc/\$p/cmdline); case \"\$c/ \$a \" in \"\$d\"/*|*\" \$d \"*|*\" \$d/\"*) exit 0;; esac; done; exit 1"; then
    echo "SKIP $r: in use on aida"; continue
  fi
  sz=$(ssh aida "du -sb ~/smol/runs/$r | cut -f1")
  t0=$(date +%s)
  ssh aida "tar -C ~/smol/runs/$(dirname $r) -cf - $(basename $r)" \
    | tee >(sha256sum | cut -d' ' -f1 > $tmp/sha) >(wc -c > $tmp/n) \
    | aws $P s3 cp - "s3://$B/$key" --storage-class GLACIER_IR --expected-size $((sz + sz / 50 + 10485760)) --only-show-errors \
    || { echo "FAIL $r (upload)"; continue; }
  sleep 1
  n=$(tr -d ' ' < $tmp/n); sha=$(cat $tmp/sha)
  got=$(aws $P s3api head-object --bucket $B --key "$key" --query ContentLength --output text)
  if [ "$n" != "$got" ] || [ -z "$sha" ]; then echo "FAIL $r: tar $n bytes, S3 $got"; continue; fi
  printf "%s  %s\n" "$sha" "$(basename $r).tar" | aws $P s3 cp - "s3://$B/manifest/$r.sha256" --only-show-errors
  dt=$(( $(date +%s) - t0 )); printf "%s\t%s\t%s\t%s\t%s\n" "$(date -u +%FT%TZ)" "$key" "$n" "$sha" "$dt" >> "$MAN"
  echo "ok $r  $n bytes  ${dt}s  $(( n / (dt > 0 ? dt : 1) / 1000000 )) MB/s"
  if [ $DEL = 1 ]; then ssh aida "rm -rf ~/smol/runs/$r" && echo "   deleted on aida"; fi
done
rm -rf $tmp
