#!/bin/bash
# Run a makedf job list locally (parallel), e.g. for small-scale tests without the grid.
# Usage: grid/run_local.sh <joblist> <outdir> [nproc]
set -e
LIST=$(readlink -f ${1:?joblist}); OUT=${2:?outdir}; NP=${3:-8}
REPO=$(cd $(dirname $(readlink -f $0))/.. && pwd)
source $REPO/makedf/scripts/setup_env.sh >/dev/null 2>&1
mkdir -p $OUT/logs
grep . $LIST | xargs -P $NP -L 1 bash -c '
  name=$0; inputs=$(echo $1 | tr "," " ")
  [ -f '"$OUT"'/$name.root ] && exit 0
  '"$REPO"'/makedf/build/makedf -c '"$REPO"'/makedf/config -o '"$OUT"'/$name.root.tmp $inputs > '"$OUT"'/logs/$name.log 2>&1 \
    && mv '"$OUT"'/$name.root.tmp '"$OUT"'/$name.root || echo "FAILED $name"'
