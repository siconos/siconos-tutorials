#-----------------------------------------------
# Ctest driver for siconos install.
# Target : continuous integration on gitlab-ci,
# aims at providing a proper install of siconos
# to test examples.
#
# Input variables :
# - SICONOS_INSTALL_DIR : where to install siconos. Default : ../install-siconos
# - USER_FILE : user option file used by cmake to configure siconos. Default : siconos_conf.cmake.
#   Warning : always searched in siconos-tutorials/ci directory.
# - OSNAME : host system name (used to qualify cdash build). If not set, try to catch info
#   using common commands (lsb_release ...)
# - WITH_TESTS : if ON, execute tests
# ----------------------------------------------

# Assumes :
# - 'siconos' git repository is available in ../
# - 'siconos-tutorials' git repository is available in ../.. (for ctest drivers)

# used as CMAKE_INSTALL_PREFIX
if(NOT SICONOS_INSTALL_DIR)
  set(SICONOS_INSTALL_DIR ../install-siconos/)
endif()

# Current testing model (Must be one of Experimental, Continuous, or Nightly)
if(NOT model)
  set(model Experimental)
endif()

if(NOT USER_FILE)
  set(USER_FILE siconos_default.cmake)
endif()

if(NOT WITH_TESTS)
  set(WITH_TESTS OFF)
endif()

# To deactivate submission to cdash. Default = submit
if(NOT DO_SUBMIT)
  set(DO_SUBMIT ON)
endif()

# For CI, we assume :
# - SICONOS_TUTORIAL_SOURCE_DIR/siconos : siconos git repo
# - SICONOS_TUTORIAL_SOURCE_DIR/examples : directory which contains examples to be tested.


# With SICONOS_TUTORIAL_SOURCE_DIR the path to git repository for siconos tutorial
get_filename_component(SICONOS_TUTORIAL_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR} DIRECTORY) # path to this (ctest_driver) file
# -- job : build and install siconos --
message("--- Start conf for siconos install.")
# - Source dir and path to siconos install
# We assume CI setup, with build dir in siconos-tutorial repository and
# siconos clone is in siconos-tutorial/siconos
# The 'cmake' source dir is siconos.
if(NOT CTEST_SOURCE_DIRECTORY)
  set(CTEST_SOURCE_DIRECTORY ../siconos)
endif()

# Build name (for cdash)
if(NOT CTEST_BUILD_NAME)
  # Get hash for commit of current version of Siconos
  file(READ ${SICONOS_TUTORIAL_SOURCE_DIR}/siconos-commit-number.txt SICO_REF)
  string(STRIP ${SICO_REF} SICO_REF)
  include(${CTEST_SOURCE_DIRECTORY}/cmake/SiconosVersion.cmake)
  set(CTEST_BUILD_NAME "Siconos install (${SICONOS_VERSION}-devel, ref commit : ${SICO_REF})")
  if(EXTRA_NAME)
    set(CTEST_BUILD_NAME "${CTEST_BUILD_NAME} - ${EXTRA_NAME}")
  endif()
endif()


set(SICONOS_CMAKE_OPTIONS -DUSER_OPTIONS_FILE=${SICONOS_TUTORIAL_SOURCE_DIR}/ci/${USER_FILE})
list(APPEND SICONOS_CMAKE_OPTIONS -DCMAKE_INSTALL_PREFIX=${SICONOS_INSTALL_DIR})
list(APPEND SICONOS_CMAKE_OPTIONS -DCMAKE_CXX_STANDARD=11)
list(APPEND SICONOS_CMAKE_OPTIONS -DSICONOS_USE_BOOST_FOR_CXX11=OFF)
list(APPEND SICONOS_CMAKE_OPTIONS -DWITH_TESTING=${WITH_TESTS})
set(current_project siconos_install)
cmake_host_system_information(RESULT NP QUERY NUMBER_OF_LOGICAL_CORES)
set(CTEST_BUILD_FLAGS -j${NP})
# Parallel build only for siconos_install. For examples it leads to: warning: jobserver unavailable: using -j1. Add `+' to parent make rule.
set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE ${CTEST_SOURCE_DIRECTORY}/cmake/valgrind.supp)

include(${SICONOS_TUTORIAL_SOURCE_DIR}/ci/ctest_common.cmake)
