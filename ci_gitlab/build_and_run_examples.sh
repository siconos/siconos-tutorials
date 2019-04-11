#!bin/bash
# This scripts configures, builds and executes tests from
# cmake config defined in siconos-tutorials/examples.
#
# It assumes a proper install of Siconos in ${CI_PROJECT_DIR}/install-siconos
# 
: ${CI_PROJECT_DIR:?"Please set environment variable CI_PROJECT_DIR with 'siconos-tutorials' repository (absolute) path."}
: ${SICONOS_INSTALL_DIR:?"Please set environment variable SICONOS_INSTALL_DIR  with Siconos installation path (absolute)."}

# -- install some extra dependencies --
pip3 install -U -r ./ci_gitlab/requirements.txt > /dev/null

# -- Create build dir --
mkdir -p $CI_PROJECT_DIR/build/examples
cd $CI_PROJECT_DIR/build/examples

# -- Run ctest --
# It will configure, build, and test siconos examples.
# 

ctest -S ${CI_PROJECT_DIR}/ci_gitlab/ctest_driver_examples.cmake -Dmodel=$ctest_build_model -DSICONOS_INSTALL_DIR=${SICONOS_INSTALL_DIR} -DOSNAME=$IMAGE_NAME -DCDASH_SUBMIT=$cdash_submit -V 

#ctest -S ${ref_path}/ci/ctest_driver_examples.cmake  -DSICONOS_INSTALL_DIR=${ref_path}/install-siconos -Dmodel=$CTEST_MODEL -DOSNAME=$1 -DEXTRA_NAME="$EXTRA_BUILDNAME" -V


