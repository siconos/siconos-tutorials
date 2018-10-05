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
if(NOT WITH_TESTS)
  # tests on for examples but off for siconos install. 
  set(WITH_TESTS ON)
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


# For CI, we assume :
# - SICONOS_TUTORIAL_SOURCE_DIR/siconos : siconos git repo
# - SICONOS_TUTORIAL_SOURCE_DIR/examples : directory which contains examples to be tested.


# With SICONOS_TUTORIAL_SOURCE_DIR the path to git repository for siconos tutorial
get_filename_component(SICONOS_TUTORIAL_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR} DIRECTORY) # path to this (ctest_driver) file

# ================= Start config which depends on job's type (siconos_install or examples) ====================


# -- job : build and install siconos --
if(JOB_NAME STREQUAL "siconos_install")
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
    set(CTEST_BUILD_NAME "Siconos install for examples")
  endif()
  set(SICONOS_CMAKE_OPTIONS -DUSER_OPTIONS_FILE=${SICONOS_TUTORIAL_SOURCE_DIR}/ci/siconos_conf.cmake)
  list(APPEND SICONOS_CMAKE_OPTIONS -DCMAKE_INSTALL_PREFIX=${SICONOS_INSTALL_DIR})
  list(APPEND SICONOS_CMAKE_OPTIONS -DCMAKE_CXX_STANDARD=11)
  list(APPEND SICONOS_CMAKE_OPTIONS -DSICONOS_USE_BOOST_FOR_CXX11=OFF)
  set(current_project siconos)
  set(CTEST_BUILD_FLAGS -j${NP})
  # Parallel build only for siconos_install. For examples it leads to ‘warning: jobserver unavailable: using -j1. Add `+' to parent make rule.’

elseif(JOB_NAME STREQUAL "examples")
  # -- job : build and test examples -- 
  message("--- Start conf for siconos examples build and tests.")
  
  # Get hash for commit of installed version of Siconos
  file(READ ${SICONOS_INSTALL_DIR}/siconos-commit-number.txt SICO_REF)
  string(STRIP ${SICO_REF} SICO_REF)
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

  set(SICONOS_CMAKE_OPTIONS -Dsiconos_DIR=${SICONOS_INSTALL_DIR}/share/siconos/cmake/)
  set(current_project siconos_examples)

endif()
# ================= End config which depends on job's type (siconos_install or examples) ====================

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



# With gitlab-ci, runner name is too long and useless ...
string(FIND ${hostname} "runner-" on_ci) 
if(on_ci GREATER -1)
  set(hostname "gitlab-ci runner on $ENV{CI_RUNNER_DESCRIPTION}")
endif()

if(JOB_NAME STREQUAL "examples")
  set(hostname "${hostname}, based on Siconos commit ${SICO_REF}")
endif()


# Host description
if(NOT OSNAME)
  set(OSNAME ${osname}) # Use -DOSNAME=docker_image name on CI
endif()
if(NOT CTEST_SITE)
  set(CTEST_SITE "${OSNAME} ${osrelease}, ${osplatform}, ${hostname}")
  #set(CTEST_SITE "${OSNAME} ${osrelease}, ${osplatform}")#, ${hostname}")
endif()

ctest_start(${model})

set(CTEST_CONFIGURE_COMMAND ${CMAKE_COMMAND})
foreach(option ${SICONOS_CMAKE_OPTIONS})
  set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${option}")
endforeach()

set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${CTEST_SOURCE_DIRECTORY}")

ctest_configure()

# --- Build ---

if(NOT CTEST_BUILD_CONFIGURATION)
  set(CTEST_BUILD_CONFIGURATION "Profiling")
endif()
ctest_build(
 PROJECT_NAME ${current_project}
 )


# -- Tests --
if(WITH_TESTS)

  ctest_test(
    PARALLEL_LEVEL NP
    RETURN_VALUE TEST_RETURN_VAL
    )
endif()

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


# ============= Summary =============
message(STATUS "\n============================================ Summary ============================================")
message(STATUS "CTest process for ${curren_project} (job name = ${JOB_NAME}) has ended.")
message(STATUS "Ctest model is: ${model}")
message(STATUS "Ctest executed on sources directory : ${CTEST_SOURCE_DIRECTORY}")
message(STATUS "Build name (cdash) : ${CTEST_BUILD_NAME}")
message(STATUS "Site (cdash) : ${CTEST_SITE}")
message(STATUS "=================================================================================================\n")

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

