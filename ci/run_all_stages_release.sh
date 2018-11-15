#
#
# 1 - clone, build and install last versions of oce and pythonocc
# 2 - clone, build and install siconos release 'tag' (input arg $2)
# 3 - build and tests examples of release 'tag'
# 
# Note : this script takes osname (from docker image in gitlab-ci script) and tag value as args.
# 
if  [ -x "$(command -v nproc)" ]; then
   export nbprocs=`nproc --all`  # linux
elif  [ -x "$(command -v sysctl)" ]; then
   export nbprocs=`sysctl -n hw.ncpu` # macos
else
   export nbprocs=2
fi
# Check if CI_PROJECT_DIR is set AND not empty
: ${CI_PROJECT_DIR:?"Please set environment variable CI_PROJECT_DIR with 'siconos-tutorials' repository (absolute) path."}
# Create directory to install libs.
pip3 install -U matplotlib
# -- Siconos clone,  build, install --
git clone https://github.com/siconos/siconos.git
cd siconos
git checkout tags/$2
git rev-parse --short HEAD > ${CI_PROJECT_DIR}/siconos-commit-number.txt
sed -i 's/bipop/tripop/g' ./CTestConfig.cmake # cdash site has changed since 4.2 ...
sed -i 's/http/https/g' ./CTestConfig.cmake
cp ./CTestConfig.cmake ./examples/
mkdir ${CI_PROJECT_DIR}/build-siconos
cd ${CI_PROJECT_DIR}/build-siconos
export buildname="Siconos install  (release/tag $tag, with OCE)"
ctest -S ${CI_PROJECT_DIR}/ci/ctest_driver_install_siconos.cmake -Dmodel=$CTEST_MODEL -DSICONOS_INSTALL_DIR=${CI_PROJECT_DIR}/install-siconos -DOSNAME=$1 -DUSER_FILE=siconos4-2-0_with_mechanisms.cmake -DCTEST_BUILD_NAME="$buildname" -V
make install -j $nbprocs > /dev/null
mv ${CI_PROJECT_DIR}/siconos-commit-number.txt ${CI_PROJECT_DIR}/install-siconos/

# -- build, tests examples --
echo "-------------> Build and run examples"
mkdir ${CI_PROJECT_DIR}/build-examples
cd ${CI_PROJECT_DIR}/build-examples
export buildname="Siconos examples (release/tag $tag, with OCE)"
export PATH=$CI_PROJECT_DIR/install-siconos/bin:$PATH
ctest -S ${CI_PROJECT_DIR}/ci/ctest_driver_examples.cmake -DSICONOS_INSTALL_DIR=${CI_PROJECT_DIR}/install-siconos -Dmodel=$CTEST_MODEL -DOSNAME=$1 -DCTEST_BUILD_NAME="$buildname" -DCTEST_SOURCE_DIRECTORY=${CI_PROJECT_DIR}/siconos/examples -V
