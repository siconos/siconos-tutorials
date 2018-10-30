#!bin/bash
# Note : this script takes osname (from docker image in gitlab-ci script) as arg.

pip3 install -U -r ./ci/requirements.txt
export ref_path=$PWD
mkdir build-examples
cd build-examples
ctest -S ${ref_path}/ci/ctest_driver_examples.cmake  -DSICONOS_INSTALL_DIR=${ref_path}/install-siconos -Dmodel=Continuous -DOSNAME=$1 -VV


