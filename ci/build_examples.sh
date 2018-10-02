#!bin/bash
pip3 install -U -r ./ci/requirements.txt
mkdir build-examples
cd build-examples
export siconos_DIR=../install-siconos/share/siconos/cmake
cmake ../examples
make -j 4
