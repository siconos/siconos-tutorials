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
mkdir oce-last pythonocc siconos
cd oce-last
cmake ../../oce -DOCE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=Release
make -j $nbprocs
make install
cd ../pythonocc
cmake ../../pythonocc-core -DCMAKE_BUILD_TYPE=Release
make install -j $nbprocs
cd $ref_path
python3 -c 'import OCC'

# h5py required for mechanics tests ...
pip3 install -U -r ./ci/requirements.txt

git clone https://github.com/siconos/siconos.git
# Get last commit id, will be used for buildname on cdash.
cd siconos
git rev-parse --short HEAD > ${ref_path}/siconos-commit-number.txt
#
cd ${ref_path}/build/siconos
ctest -S ${ref_path}/ci/ctest_driver.cmake -V -DJOB_NAME=siconos_install -Dmodel=Continuous -DSICONOS_INSTALL_DIR=${ref_path}/install-siconos -DOSNAME=$1 -DUSER_FILE=siconos_with_mechanisms.cmake -DCTEST_BUILD_NAME="Siconos (with mechanisms) install for examples" -DCTEST_SOURCE_DIRECTORY=${ref_path}/siconos

#make -j ${nbprocs}
make install
mv ${ref_path}/siconos-commit-number.txt ${ref_path}/install-siconos/
