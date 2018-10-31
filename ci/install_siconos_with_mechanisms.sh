#!bin/bash

# Full install of siconos (all components, with oce).

# 1 - clone, build and install last versions of oce and pythonocc
# 2 - clone, build and install siconos

# Note : this script takes osname (from docker image in gitlab-ci script) as arg.


# Get number of procs
if  [ -x "$(command -v nproc)" ]; then
   export nbprocs=`nproc --all`  # linux
elif  [ -x "$(command -v sysctl)" ]; then
   export nbprocs=`sysctl -n hw.ncpu` # macos
else
   export nbprocs=2
fi

export ref_path=$PWD
env

# -- occ/python occ install --
mkdir build
cd build
sh ${ref_path}/ci/install_oce.sh $2
cd $ref_path

# -- siconos install --
git clone https://github.com/siconos/siconos.git
# Get last commit id, will be used for buildname on cdash.
cd siconos
git rev-parse --short HEAD > ${ref_path}/siconos-commit-number.txt
#
cd ${ref_path}
mkdir build-siconos
cd build-siconos
ctest -S ${ref_path}/ci/ctest_driver_install_siconos.cmake -V -Dmodel=Continuous -DSICONOS_INSTALL_DIR=${ref_path}/install-siconos -DOSNAME=$1 -DUSER_FILE=siconos_with_mechanisms.cmake -DEXTRA_NAME="With mechanisms/OCE" -DWITH_TESTS=ON -DCTEST_BUILD_CONFIGURATION=Release

#make -j ${nbprocs}
make install
mv ${ref_path}/siconos-commit-number.txt ${ref_path}/install-siconos/
