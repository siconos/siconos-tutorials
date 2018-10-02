#!bin/bash
pip3 install -U -r ./ci/requirements.txt
git clone https://github.com/siconos/siconos.git
mkdir build
cd build
cmake ../siconos -DUSER_OPTIONS_FILE=$PWD/../ci/siconos_conf.cmake -DCMAKE_INSTALL_PREFIX=../install-siconos
make -j 4
make install
