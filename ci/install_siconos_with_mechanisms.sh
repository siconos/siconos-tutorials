#!bin/bash

# Full install of siconos (all components, with oce).

# 1 - clone, build and install last versions of oce and pythonocc
# 2 - clone, build and install siconos

# Note : this script takes osname (from docker image in gitlab-ci script) as arg.

# Check if CI_PROJECT_DIR is set AND not empty
: ${CI_PROJECT_DIR:?"Please set environment variable CI_PROJECT_DIR with 'siconos-tutorials' repository (absolute) path."}
# -- siconos download, build, test and install --
git clone https://github.com/siconos/siconos.git
# Get last commit id, will be used for buildname on cdash.
cd siconos
git rev-parse --short HEAD > ${CI_PROJECT_DIR}/siconos-commit-number.txt
#
cd ${CI_PROJECT_DIR}
mkdir -p build-siconos
cd build-siconos
ctest -S ${CI_PROJECT_DIR}/ci/ctest_driver_install_siconos.cmake -Dmodel=$CTEST_MODEL -DSICONOS_INSTALL_DIR=${CI_PROJECT_DIR}/install-siconos -DOSNAME=$1 -DUSER_FILE=siconos_with_mechanisms.cmake -DEXTRA_NAME="With mechanisms/OCE" -V

make install > /dev/null
mv ${CI_PROJECT_DIR}/siconos-commit-number.txt ${CI_PROJECT_DIR}/install-siconos/
