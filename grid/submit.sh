#!/bin/bash
# Submit makedf grid jobs.
# Usage: grid/submit.sh <tag> <joblist> [extra jobsub args]
#   tag      : DF production tag (tarball grid/tarballs/makedf_<tag>.tar.gz must exist)
#   joblist  : file from make_joblists.py (one job per line)
# Outputs go to /pnfs/sbnd/scratch/users/$USER/MvA-Loop/df/<tag>/<joblist-basename>/
set -e
TAG=${1:?usage: submit.sh <tag> <joblist>}
LIST=$(readlink -f ${2:?usage: submit.sh <tag> <joblist>})
shift 2
REPO=$(cd $(dirname $(readlink -f $0))/.. && pwd)
TAR=$REPO/grid/tarballs/makedf_${TAG}.tar.gz
[ -f $TAR ] || { echo "missing $TAR; run grid/build_tarball.sh $TAG"; exit 1; }
NAME=${OUTNAME:-$(basename $LIST .txt)}  # OUTNAME overrides (e.g. for resubmissions)
OUTDIR=/pnfs/sbnd/scratch/users/${USER}/MvA-Loop/df/${TAG}/${NAME}
mkdir -p $OUTDIR/logs
N=$(grep -c . $LIST)
echo "submitting $N jobs: $LIST -> $OUTDIR"
jobsub_submit -G sbnd -N $N \
  --memory=2500MB --disk=8GB --expected-lifetime=8h \
  --resource-provides=usage_model=${USAGE_MODEL:-OPPORTUNISTIC,DEDICATED,OFFSITE} \
  --singularity-image /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-el9:latest \
  --tar_file_name dropbox://$TAR \
  -f dropbox://$LIST \
  -e JOBLIST=$(basename $LIST) -e OUTDIR=$OUTDIR -e TAG=$TAG \
  "$@" \
  file://$REPO/grid/gridjob.sh
