# Environment for building/running makedf (and PyROOT aux converters).
# Same spack ROOT on interactive nodes and on the grid (el9 container + cvmfs).
source /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh
spack load root@6.28.12 arch=linux-almalinux9-x86_64_v3
spack load gcc@12.2.0 arch=linux-almalinux9-x86_64_v3
spack load cmake@3.27.9%gcc@11.4.1 arch=linux-almalinux9-x86_64_v3
export LD_LIBRARY_PATH=${ROOTSYS}/lib/root:${LD_LIBRARY_PATH}
