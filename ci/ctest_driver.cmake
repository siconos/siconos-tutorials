#-----------------------------------------------
# Ctest driver for siconos install.
# Target : continuous integration on gitlab-ci,
# aims at providing a proper install of siconos
# to test examples.
# ----------------------------------------------


# JOB_NAME:
# - siconos_install to build and install Siconos.
#   Assumes 'siconos' git repository is available in ../
# - examples to build and test examples of the present repository.
#   Assumes a proper install of Siconos.

if(NOT JOB_NAME)
  set(JOB_NAME "examples")
endif()

# Path to siconos install
#  if JOB_NAME is siconos_install, used as CMAKE_INSTALL_PREFIX
#  and if JOB_NAME is examples, used to find siconos.
if(NOT SICONOS_INSTALL_DIR)
  set(SICONOS_INSTALL_DIR ../install-siconos/)
endif()

# Current testing model (Must be one of Experimental, Continuous, or Nightly)
if(NOT model)
  set(model Experimental)
endif()



if(JOB_NAME STREQUAL "siconos_install")
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
    set(CTEST_BUILD_NAME "Siconos install for examples")
  endif()
  set(CTEST_CONFIG_OPTIONS "-DUSER_OPTIONS_FILE=$PWD/../ci/siconos_conf.cmake -DCMAKE_INSTALL_PREFIX=${SICONOS_INSTALL_DIR}")
  set(current_project siconos)
elseif(JOB_NAME STREQUAL "examples")
  message("--- Start conf for siconos examples build and tests.")

  # Get hash for commit of installed version of Siconos
  file(READ ${SICONOS_INSTALL_DIR}/siconos-commit-number.txt SICO_REF)
  # - Source dir and path to siconos install
  # We assume CI setup, with build dir in siconos-tutorial repository and
  # siconos install in siconos-tutorial/install-siconos.
  # The 'cmake' source dir is examples.
  if(NOT CTEST_SOURCE_DIRECTORY)
    set(CTEST_SOURCE_DIRECTORY ../examples)
  endif()
  # Build name (for cdash)
  if(NOT CTEST_BUILD_NAME)
    set(CTEST_BUILD_NAME "Siconos examples")
  endif()

  set(CTEST_CONFIG_OPTIONS "-Dsiconos_DIR=${SICONOS_INSTALL_DIR}/share/siconos/cmake/")
  set(current_project siconos_examples)

endif()

# --- Configure setup ---

# - Top level build directory -
# If not specified : current dir.
if(NOT CTEST_BINARY_DIRECTORY)
  set(CTEST_BINARY_DIRECTORY .)
endif()


if(NOT CTEST_CMAKE_GENERATOR)
  set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
endif()

# -- Query host system information --
#include(cmake_host_system_information)
cmake_host_system_information(RESULT NP QUERY NUMBER_OF_LOGICAL_CORES)
cmake_host_system_information(RESULT hostname QUERY HOSTNAME)
cmake_host_system_information(RESULT fqdn QUERY FQDN)

if(${CMAKE_VERSION} VERSION_GREATER "3.10.3") 
  cmake_host_system_information(RESULT osname QUERY OS_NAME)
  cmake_host_system_information(RESULT osrelease QUERY OS_RELEASE)
  cmake_host_system_information(RESULT osversion QUERY OS_VERSION)
  cmake_host_system_information(RESULT osplatform QUERY OS_PLATFORM)
else()
  set(osname ${CMAKE_SYSTEM_NAME})
  set(osversion ${CMAKE_SYSTEM_VERSION})
  set(osplatform ${CMAKE_SYSTEM_PROCESSOR})
endif()



# Runner name is too long and useless ...
# if(ENV{CI_RUNNER_DESCRIPTION})
#   # If on a gitlab-ci runner ...
#   set(hostname "gitlab-ci runner on $ENV{CI_RUNNER_DESCRIPTION}")
# endif()
string(FIND ${hostname} "runner-" on_ci) 
if(on_ci GREATER -1)
  set(hostname "gitlab-ci runner on $ENV{CI_RUNNER_DESCRIPTION}\
                based on Siconos commit ${SICO_REF}.")
endif()

# Host description
if(NOT CTEST_SITE)
  set(CTEST_SITE "${hostname}, ${osname}, ${osrelease}, ${osplatform}")
endif()

ctest_start(${model})

ctest_configure(OPTIONS ${CTEST_CONFIG_OPTIONS})

# --- Build ---

if(NOT CTEST_BUILD_CONFIGURATION)
  set(CTEST_BUILD_CONFIGURATION "Profiling")
endif()

ctest_build(
 PROJECT_NAME ${current_project}
 )

# -- Tests --

# check number of cores available on the host
ctest_test(
 PARALLEL_LEVEL NP
 RETURN_VALUE TEST_RETURN_VAL
 )

# -- memory check --
if(CTEST_BUILD_CONFIGURATION MATCHES "Profiling")
  find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)
  set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--quiet --leak-check=full --show-leak-kinds=definite,possible --track-origins=yes --error-limit=no --gen-suppressions=all") 
  set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--quiet --leak-check=full --show-reachable=yes --error-limit=no --gen-suppressions=all") 
  set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE ${CTEST_SOURCE_DIRECTORY}/cmake/valgrind.supp)
  
  ctest_memcheck(PARALLEL_LEVEL NP)
endif()

# -- Submission to cdash --
#ctest_submit(RETURN_VALUE SUBMIT_RETURN_VAL)

# submit failed? 
if(NOT SUBMIT_RETURN_VAL EQUAL 0)
  message(WARNING " *** submission failure *** ")
endif()

# tests failed? 
if(NOT TEST_RETURN_VAL EQUAL 0)
  message(FATAL_ERROR " *** test failure *** ")
endif()


# message(STATUS "Siconos CTest driver")
# message(STATUS "MODE is: ${MODE}")
# message(STATUS "CTEST_SOURCE_DIRECTORY is: ${CTEST_SOURCE_DIRECTORY}")
# message(STATUS "CTEST_BINARY_DIRECTORY is: ${CTEST_BINARY_DIRECTORY}")
# message(STATUS "CTEST_MODULE_PATH is: ${CTEST_MODULE_PATH}")
# message(STATUS "CTEST_BUILD_CONFIGURATION is: ${CTEST_BUILD_CONFIGURATION}")
# #######################################################################
# # this usually fails for some reasons and ctest may returns a fail code.
# # ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY}/)
# # cf discussions here:
# # https://gitlab.kitware.com/cmake/cmake/issues/17000

# if(CTEST_BINARY_DIRECTORY)
#   file(REMOVE_RECURSE ${CTEST_BINARY_DIRECTORY})
# endif()

#find_program(CTEST_GIT_COMMAND NAMES git)
#find_program(CTEST_COVERAGE_COMMAND NAMES gcov)

# if(TEST_TIMEOUT)
#   set(CTEST_TEST_TIMEOUT ${TEST_TIMEOUT})
# endif(TEST_TIMEOUT)

