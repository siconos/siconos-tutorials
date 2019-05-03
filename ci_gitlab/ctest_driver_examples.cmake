#-----------------------------------------------
# Ctest driver to build and tests examples.
# Target : continuous integration on gitlab-ci,
#
# Input variables :
# - SICONOS_INSTALL_DIR : path to Siconos install.
# - OSNAME : host system name (used to qualify cdash build). If not set, try to catch info
#   using common commands (lsb_release ...)
# ----------------------------------------------

# Assumes :
# - a proper install of siconos in SICONOS_INSTALL_DIR
# - an environment variable CI_PROJECT_DIR which contains the path
# - a file siconos-commit-number.txt in $HOME.
#   to siconos-tutorials repository.


# -- job : build and test examples -- 
message("--- Start conf for siconos examples build and tests.")

# - Source dir and path to siconos-tutorials/examples
if(NOT CTEST_SOURCE_DIRECTORY)
  set(CTEST_SOURCE_DIRECTORY $ENV{CI_PROJECT_DIR}/examples)
endif()

# Build name (for cdash)
if(NOT CTEST_BUILD_NAME)
  # Get hash for commit of installed version of Siconos
  file(READ $ENV{HOME}/siconos-commit-number.txt SICO_REF)
  string(STRIP ${SICO_REF} SICO_REF)
  set(CTEST_BUILD_NAME "Examples (based on Siconos commit : ${SICO_REF})")
  if(EXTRA_NAME)
    set(CTEST_BUILD_NAME "${CTEST_BUILD_NAME} - ${EXTRA_NAME}.")
  endif()
endif()

include($ENV{CI_PROJECT_DIR}/examples/cmake/SiconosRequiredVersion.cmake)
set(ConfigPackageLocation lib/cmake/siconos-${SICONOS_REQUIRED_VERSION})
set(SICONOS_CMAKE_OPTIONS -Dsiconos_DIR=${SICONOS_INSTALL_DIR}/${ConfigPackageLocation})
set(current_project siconos_examples)
set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE ${SICONOS_INSTALL_DIR}/${ConfigPackageLocation})
include($ENV{CI_PROJECT_DIR}/ci_gitlab/ctest_common.cmake)


