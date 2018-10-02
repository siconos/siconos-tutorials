#!bin/bash
apt update -qq && apt install -y -qq cmake git-core wget make \
			      libboost-dev libgmp-dev swig gcc gfortran g++ liblapack-dev libatlas-base-dev \
			      lp-solve liblpsolve55-dev python3-dev libpython3-dev bash swig doxygen python3-dev python3-pip graphviz htop

pip3 install -U -r ./ci/requirements.txt
git clone https://github.com/siconos/siconos.git
mkdir build
cd build
cmake ../siconos -DUSER_OPTIONS_FILE=$PWD/../ci/siconos_conf.cmake -DCMAKE_INSTALL_PREFIX=../install-siconos
make -j 4
make install
