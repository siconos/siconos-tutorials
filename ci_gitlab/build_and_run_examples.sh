#!bin/bash
# This scripts configures, builds and executes tests from
# cmake config defined in siconos-tutorials/examples.
#
# It assumes a proper install of Siconos in ${CI_PROJECT_DIR}/install-siconos
# 
: ${CI_PROJECT_DIR:?"Please set environment variable CI_PROJECT_DIR with 'siconos-tutorials' repository (absolute) path."}

# Default build dir, if not set
BUILD_DIR="${BUILD_DIR:=${CI_PROJECT_DIR}/build-examples}" 
# Default ctest mode
CTEST_BUILD_MODEL="${CTEST_BUILD_MODEL:=Experimental}"
# Set to 1 to allow -jN, 0 to restrict to -j1.
PARALLEL_BUILD="${PARALLEL_BUILD=:=1}"
# Default: submit to cdash
CDASH_SUBMIT="${CDASH_SUBMIT=:=1}"


# -- Run ctest --
# It will configure, build, and test siconos examples.
ctest -S ${CI_PROJECT_DIR}/ci_gitlab/ctest_driver_examples.cmake -Dmodel=$CTEST_BUILD_MODEL -DALLOW_PARALLEL_BUILD=$PARALLEL_BUILD -DCDASH_SUBMIT=$CDASH_SUBMIT -V --output-junit test_results.xml -DCTEST_BINARY_DIRECTORY=${BUILD_DIR} -DCTEST_SOURCE_DIRECTORY=$CI_PROJECT_DIR/examples


