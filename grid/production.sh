#!/bin/bash
# Build the tarball and job lists for a DF production tag and submit everything.
# Usage: grid/production.sh <tag> [samples...]
#   samples default: me_fhc_mc me_fhc_data   (also: me_rhc_mc me_rhc_data)
# Small-scale tests: set PLAYLISTS=1A MAXFILES=10 in the environment.
set -e
TAG=${1:?usage: production.sh <tag> [samples]}
shift
SAMPLES=${@:-me_fhc_mc me_fhc_data}
REPO=$(cd $(dirname $(readlink -f $0))/.. && pwd)
cd $REPO
grid/build_tarball.sh $TAG
mkdir -p joblists/$TAG
for s in $SAMPLES; do
  case $s in
    me_fhc_mc)   args="--sample MediumEnergy_FHC --kind MC" ;;
    me_fhc_data) args="--sample MediumEnergy_FHC --kind Data" ;;
    me_rhc_mc)   args="--sample MediumEnergy_RHC --kind MC" ;;
    me_rhc_data) args="--sample MediumEnergy_RHC --kind Data" ;;
    le_fhc_mc)   args="--sample LowEnergy_FHC --kind MC" ;;
    le_fhc_data) args="--sample LowEnergy_FHC --kind Data" ;;
    *) echo "unknown sample $s"; exit 1 ;;
  esac
  [ -n "$PLAYLISTS" ] && args="$args --playlists $PLAYLISTS"
  [ -n "$MAXFILES" ] && [[ $s == *_mc ]] && args="$args --max-files $MAXFILES"
  python3 grid/make_joblists.py $args -o joblists/$TAG/$s.txt
  grid/submit.sh $TAG joblists/$TAG/$s.txt
done
