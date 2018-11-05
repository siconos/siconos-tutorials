#!bin/bash

# Standard install of siconos (all components, without oce).

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
ctest -S ${ref_path}/ci/ctest_driver_install_siconos.cmake -Dmodel=Continuous -DSICONOS_INSTALL_DIR=${ref_path}/install-siconos -DOSNAME=$1

make install
mv ${ref_path}/siconos-commit-number.txt ${ref_path}/install-siconos/
