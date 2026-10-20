#!/bin/bash
# Build makedf and package binary + configs + job script for the grid.
# Usage: grid/build_tarball.sh <tag>      -> grid/tarballs/makedf_<tag>.tar.gz
set -e
TAG=${1:?usage: build_tarball.sh <tag>}
REPO=$(cd $(dirname $(readlink -f $0))/.. && pwd)
source $REPO/makedf/scripts/setup_env.sh >/dev/null 2>&1
cmake -S $REPO/makedf -B $REPO/makedf/build -DCMAKE_CXX_COMPILER=$(which g++) >/dev/null
cmake --build $REPO/makedf/build -j8
STAGE=$(mktemp -d)
mkdir -p $STAGE/makedf/bin $STAGE/makedf/scripts
cp $REPO/makedf/build/makedf $STAGE/makedf/bin/
cp -r $REPO/makedf/config $STAGE/makedf/
cp $REPO/makedf/scripts/setup_env.sh $STAGE/makedf/scripts/
cp $REPO/grid/gridjob.sh $STAGE/
# aux inputs for the Geant hadron reweight (MAT data/hadronReweight + per-playlist renormalization)
AUX=${MAKEDF_AUX:-/exp/sbnd/data/users/gputnam/MvA-Loop/aux}/MATFluxAndReweightFiles
mkdir -p $STAGE/makedf/aux/MATFluxAndReweightFiles/hadronReweight $STAGE/makedf/aux/MATFluxAndReweightFiles/mhrwKineRenorm
cp $AUX/hadronReweight/*.root $STAGE/makedf/aux/MATFluxAndReweightFiles/hadronReweight/
cp $AUX/mhrwKineRenorm/Renorm_Kine_Truth_Minerva*_v10_HD_withElastics_RwgtNeutronCV_z5991-8408_a850.root \
   $STAGE/makedf/aux/MATFluxAndReweightFiles/mhrwKineRenorm/
(cd $REPO && git rev-parse HEAD 2>/dev/null; git -C $REPO status --porcelain | head -20) > $STAGE/makedf/VERSION
mkdir -p $REPO/grid/tarballs
tar czf $REPO/grid/tarballs/makedf_${TAG}.tar.gz -C $STAGE .
rm -rf $STAGE
echo "built $REPO/grid/tarballs/makedf_${TAG}.tar.gz"
