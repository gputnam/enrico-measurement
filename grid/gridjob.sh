#!/bin/bash
# Grid worker script: runs makedf on line ($PROCESS+1) of the job list.
# Environment (set by submit.sh): JOBLIST (basename), OUTDIR, TAG
LOG="makedf_${CLUSTER}_${PROCESS}.log"
exec > >(tee -a $LOG) 2>&1
echo "host $(hostname) cluster ${CLUSTER} process ${PROCESS} tag ${TAG}"
date

TARDIR=${INPUT_TAR_DIR_LOCAL:-$PWD}
LIST=${CONDOR_DIR_INPUT}/${JOBLIST}
LINE=$(sed -n "$((PROCESS + 1))p" $LIST)
OUTNAME=$(echo $LINE | awk '{print $1}')
INPUTS=$(echo $LINE | awk '{print $2}' | tr ',' ' ')
if [ -z "$OUTNAME" ]; then echo "no job line for process $PROCESS"; exit 1; fi
echo "outname $OUTNAME"
echo "inputs $INPUTS"

export MAKEDF_AUX=$TARDIR/makedf/aux
source $TARDIR/makedf/scripts/setup_env.sh
spack load ifdhc@2.8.0 arch=linux-almalinux9-x86_64_v3
export IFDH_CP_MAXRETRIES=3
cat $TARDIR/makedf/VERSION

OUT=${OUTNAME}.root
rc=1
for attempt in 1 2 3; do
  /usr/bin/time -v $TARDIR/makedf/bin/makedf -c $TARDIR/makedf/config -t "$TAG" -o $OUT $INPUTS
  rc=$?
  [ $rc -eq 0 ] && break
  echo "makedf failed rc=$rc (attempt $attempt)"; sleep 60
done
date

ifdh mkdir_p $OUTDIR/logs
if [ $rc -eq 0 ]; then
  ifdh cp $OUT $OUTDIR/$OUT
  rc=$?
fi
ifdh cp $LOG $OUTDIR/logs/${OUTNAME}.log
exit $rc
