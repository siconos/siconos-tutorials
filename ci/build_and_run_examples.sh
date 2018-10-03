#!bin/bash
pip3 install -U -r ./ci/requirements.txt
mkdir build-examples
cd build-examples
ctest -S ../ci/ctest_driver.cmake  -V -DSICONOS_INSTALL_DIR=../install-siconos -Dmodel=Continuous


