#!bin/bash
# Note : this script takes osname (from docker image in gitlab-ci script) as arg.

export ref_path=$PWD
git clone https://github.com/siconos/siconos.git
# Get last commit id, will be used for buildname on cdash.
cd siconos
git rev-parse --short HEAD > ${ref_path}/siconos-commit-number.txt
#
cd ..
mkdir build
cd build
ctest -S ${ref_path}/ci/ctest_driver.cmake -V -DJOB_NAME=siconos_install -Dmodel=Continuous -DWITH_TESTS=OFF -DSICONOS_INSTALL_DIR=${ref_path}/install-siconos -DOSNAME=$1

#cmake ../siconos -DUSER_OPTIONS_FILE=$PWD/../ci/siconos_conf.cmake -DCMAKE_INSTALL_PREFIX=../install-siconos
#if  [ -x "$(command -v nproc)" ]; then
#    export nbprocs=`nproc --all`  # linux
#elif  [ -x "$(command -v sysctl)" ]; then
#    export nbprocs=`sysctl -n hw.ncpu` # macos
#else
#    export nbprocs=2
#fi

#make -j ${nbprocs}
make install
mv ${ref_path}/siconos-commit-number.txt ${ref_path}/install-siconos/
