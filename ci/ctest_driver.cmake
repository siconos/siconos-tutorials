# a ctest driver for Experimental,Continuous,Nightly builds

# Path to siconos install 
if(NOT SICONOS_INSTALL_DIR)
  set(SICONOS_INSTALL_DIR ../install-siconos/)
endif()


# Current testing model (Must be one of Experimental, Continuous, or Nightly)
if(NOT model)
  set(model Experimental)
endif()

# --- Configure setup ---

# - Top level build directory -
# If not specified : current dir.
if(NOT CTEST_BINARY_DIRECTORY)
  set(CTEST_BINARY_DIRECTORY .)
endif()

# - Source dir and path to siconos install
# We assume CI setup, with build dir in siconos-tutorial repository and
# siconos install in siconos-tutorial/install-siconos.
# The 'cmake' source dir is examples.
if(NOT CTEST_SOURCE_DIRECTORY)
  set(CTEST_SOURCE_DIRECTORY ../examples)
endif()
# Path to siconos install 
if(NOT SICONOS_INSTALL_DIR)
  set(SICONOS_INSTALL_DIR ../install-siconos/)
endif()


if(NOT CTEST_CMAKE_GENERATOR)
  set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
endif()

# -- Query host system information --
#include(cmake_host_system_information)
cmake_host_system_information(RESULT NP QUERY NUMBER_OF_LOGICAL_CORES)
cmake_host_system_information(RESULT hostname QUERY HOSTNAME)
cmake_host_system_information(RESULT osname QUERY OS_NAME)
cmake_host_system_information(RESULT osrelease QUERY OS_RELEASE)
cmake_host_system_information(RESULT osversion QUERY OS_VERSION)
cmake_host_system_information(RESULT osplatform QUERY OS_PLATFORM)



message("${NP} ${hostname} ${osname} ${osrelease} ${osversion} ${osplatform}")

# Build name (for cdash)
if(NOT CTEST_BUILD_NAME)
  set(CTEST_BUILD_NAME "Siconos examples")
endif()

# Host description
if(NOT CTEST_SITE)
  set(CTEST_SITE "${hostname}, ${osname}, ${osrelease}, ${osplatform}")
endif()

ctest_start(${model})

ctest_configure(OPTIONS -Dsiconos_DIR=${SICONOS_INSTALL_DIR}/share/siconos/cmake/)

# --- Build ---

if(NOT CTEST_BUILD_CONFIGURATION)
  set(CTEST_BUILD_CONFIGURATION "Profiling")
endif()

ctest_build(
 PROJECT_NAME siconos_examples
 )

# -- Tests --

# check number of cores available on the host
#include(ProcessorCount)
#ProcessorCount(NP)
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
ctest_submit(RETURN_VALUE SUBMIT_RETURN_VAL)

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

