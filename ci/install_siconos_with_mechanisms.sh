#!bin/bash

# Full install of siconos (all components, with oce).

# 1 - clone, build and install last versions of oce and pythonocc
# 2 - clone, build and install siconos

# Note : this script takes osname (from docker image in gitlab-ci script) as arg.

export ref_path=$PWD
env
# -- siconos download, build, test and install --
git clone https://github.com/siconos/siconos.git
# Get last commit id, will be used for buildname on cdash.
cd siconos
git rev-parse --short HEAD > ${ref_path}/siconos-commit-number.txt
#
cd ${ref_path}
mkdir build-siconos
cd build-siconos
ctest -S ${ref_path}/ci/ctest_driver_install_siconos.cmake -Dmodel=Continuous -DSICONOS_INSTALL_DIR=${ref_path}/install-siconos -DOSNAME=$1 -DUSER_FILE=siconos_with_mechanisms.cmake -DEXTRA_NAME="With mechanisms/OCE"

make install
mv ${ref_path}/siconos-commit-number.txt ${ref_path}/install-siconos/
