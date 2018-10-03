#!bin/bash
#pip3 install -U scipy
#pip3 install -U pytest
env

git clone https://github.com/siconos/siconos.git
# Get last commit id, will be used for buildname on cdash.
git rev-parse --short HEAD > ./siconos-commit-number.txt
#
mkdir build
cd build
ctest -S ../ci/ctest_driver_siconos_install.cmake -V -DSICONOS_INSTALL_DIR=../install-siconos -DJOB_NAME=siconos_install

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
mv ../siconos-commit-number.txt ../install-siconos/
