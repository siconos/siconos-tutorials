#----------------------------------------------------------------
# Ctest driver : command lines common to all
# ctest drivers.
#
# Usage, add at the end of your driver file :
# include(${SICONOS_TUTORIAL_SOURCE_DIR}/ci/ctest_common.cmake)
# See for instance ctest_driver_install_siconos.cmake
# ---------------------------------------------------------------



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

# Host description
if(NOT OSNAME)
  set(OSNAME ${osname}) # Use -DOSNAME=docker_image name on CI
endif()
if(NOT CTEST_SITE)
  set(CTEST_SITE "${OSNAME} ${osrelease}, ${osplatform}, ${hostname}")
endif()

ctest_start(${model})

set(CTEST_CONFIGURE_COMMAND ${CMAKE_COMMAND})
foreach(option ${SICONOS_CMAKE_OPTIONS})
  set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${option}")
endforeach()

set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${CTEST_SOURCE_DIRECTORY}")

ctest_configure(CAPTURE_CMAKE_ERROR CONFIGURE_STATUS)

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
  ctest_memcheck(PARALLEL_LEVEL NP)
endif()

# -- Submission to cdash --
if(DO_SUBMIT)
  ctest_submit(RETURN_VALUE SUBMIT_RETURN_VAL)
  # submit failed?
  if(NOT SUBMIT_RETURN_VAL EQUAL 0)
    message(WARNING " *** submission failure *** ")
  endif()
endif()

# tests failed?
if(WITH_TESTS)
  if(NOT TEST_RETURN_VAL EQUAL 0)
    message(FATAL_ERROR " *** test failure *** ")
  endif()
endif()

# ============= Summary =============
message(STATUS "\n============================================ Summary ============================================")
message(STATUS "CTest process for ${current_project} has ended.")
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

