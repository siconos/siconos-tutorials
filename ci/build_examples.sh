mkdir build-examples
cd build-examples
export siconos_DIR=../install-siconos/share/siconos/cmake
cmake ../siconos-tutorial/examples
make -j 4
