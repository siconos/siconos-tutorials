#-----------------------------------------------
# Ctest driver to build and tests examples.
# Target : continuous integration on gitlab-ci,
#
# Input variables :
# - SICONOS_INSTALL_DIR : where to install siconos. Default : ../install-siconos
# - OSNAME : host system name (used to qualify cdash build). If not set, try to catch info
#   using common commands (lsb_release ...)
# ----------------------------------------------

# Assumes :
# - 'siconos-tutorials' git repository is available in ../.. (for ctest drivers)
# - a proper install of siconos in SICONOS_INSTALL_DIR

# Path to siconos install
if(NOT SICONOS_INSTALL_DIR)
  set(SICONOS_INSTALL_DIR ../install-siconos/)
endif()

# Current testing model (Must be one of Experimental, Continuous, or Nightly)
if(NOT model)
  set(model Experimental)
endif()

# For CI, we assume :
# - SICONOS_TUTORIAL_SOURCE_DIR/siconos : siconos git repo
# - SICONOS_TUTORIAL_SOURCE_DIR/examples : directory which contains examples to be tested.


# With SICONOS_TUTORIAL_SOURCE_DIR the path to git repository for siconos tutorial
get_filename_component(SICONOS_TUTORIAL_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR} DIRECTORY) # path to this (ctest_driver) file
# -- job : build and test examples -- 
message("--- Start conf for siconos examples build and tests.")

# - Source dir and path to siconos install
# We assume CI setup, with build dir in siconos-tutorial repository and
# siconos install in siconos-tutorial/install-siconos.
# The 'cmake' source dir is examples.
if(NOT CTEST_SOURCE_DIRECTORY)
  set(CTEST_SOURCE_DIRECTORY ../examples)
endif()

# Build name (for cdash)
if(NOT CTEST_BUILD_NAME)
  # Get hash for commit of installed version of Siconos
  file(READ ${SICONOS_INSTALL_DIR}/siconos-commit-number.txt SICO_REF)
  string(STRIP ${SICO_REF} SICO_REF)
  set(CTEST_BUILD_NAME "Examples (based on Siconos commit : ${SICO_REF})")
  if(EXTRA_NAME)
    set(CTEST_BUILD_NAME "${CTEST_BUILD_NAME} - ${EXTRA_NAME}.")
  endif()
endif()

set(SICONOS_CMAKE_OPTIONS -Dsiconos_DIR=${SICONOS_INSTALL_DIR}/share/siconos/cmake/)
set(current_project siconos_examples)
set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE ${siconos_ROOT_DIR}/share/siconos/cmake/valgrind.supp)

include(${SICONOS_TUTORIAL_SOURCE_DIR}/ci/ctest_common.cmake)
