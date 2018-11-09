#
#
# 1 - clone, build and install last versions of oce and pythonocc
# 2 - clone, build and install siconos release 'tag' (input arg $2)
# 3 - build and tests examples of release 'tag'
# 
# Note : this script takes osname (from docker image in gitlab-ci script) and tag value as args.
# 
export ref_path=$PWD
pip3 install -U matplotlib
# -- Siconos clone,  build, install --
git clone https://github.com/siconos/siconos.git
cd siconos
git checkout tags/$2
git rev-parse --short HEAD > ${ref_path}/siconos-commit-number.txt
sed -i 's/bipop/tripop/g' ./CTestConfig.cmake # cdash site has changed since 4.2 ...
mkdir ${ref_path}/build-siconos
cd ${ref_path}/build-siconos
export buildname="Siconos install (release/tag $tag, with OCE)"
ctest -S ${ref_path}/ci/ctest_driver_install_siconos.cmake -Dmodel=Continuous -DSICONOS_INSTALL_DIR=${ref_path}/install-siconos -DOSNAME=$1 -DUSER_FILE=siconos4-2-0_with_mechanisms.cmake -DCTEST_BUILD_NAME="$buildname" -V
make install
mv ${ref_path}/siconos-commit-number.txt ${ref_path}/install-siconos/

# -- build, tests examples --
mkdir ${ref_path}/build-examples
cd ${ref_path}/build-examples
export buildname="Siconos (release $tag; with OCE), run all examples"
export PATH=$ref_path/install-siconos/bin:$PATH
ctest -S ${ref_path}/ci/ctest_driver_examples.cmake -DSICONOS_INSTALL_DIR=${ref_path}/install-siconos -Dmodel=Continuous -DOSNAME=$1 -DCTEST_BUILD_NAME="$buildname" -DCTEST_SOURCE_DIRECTORY=${ref_path}/siconos/examples -V
