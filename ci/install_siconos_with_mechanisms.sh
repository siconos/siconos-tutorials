#!bin/bash
# Note : this script takes osname (from docker image in gitlab-ci script) as arg.

if  [ -x "$(command -v nproc)" ]; then
   export nbprocs=`nproc --all`  # linux
elif  [ -x "$(command -v sysctl)" ]; then
   export nbprocs=`sysctl -n hw.ncpu` # macos
else
   export nbprocs=2
fi

export ref_path=$PWD
env



# -- occ/python occ install --
git clone https://github.com/tpaviot/oce.git
git clone https://github.com/tpaviot/pythonocc-core.git
mkdir build
cd build
mkdir oce-last pythonocc
cd oce-last
cmake ../../oce -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=Release
make -j $nbprocs
make install
cd ../pythonocc
cmake ../../pythonocc-core -DCMAKE_BUILD_TYPE=Release
make install -j $nbprocs
cd $ref_path
python -c 'import OCC'

git clone https://github.com/siconos/siconos.git
# Get last commit id, will be used for buildname on cdash.
cd siconos
git rev-parse --short HEAD > ${ref_path}/siconos-commit-number.txt
#
cd ..
mkdir build
cd build
ctest -S ${ref_path}/ci/ctest_driver.cmake -V -DJOB_NAME=siconos_install -Dmodel=Continuous -DWITH_TESTS=ON -DSICONOS_INSTALL_DIR=${ref_path}/install-siconos -DOSNAME=$1 -DUSER_FILE=siconos_with_mechanisms.cmake -DCTEST_BUILD_NAME="Siconos with mechanisms"

#make -j ${nbprocs}
make install
mv ${ref_path}/siconos-commit-number.txt ${ref_path}/install-siconos/