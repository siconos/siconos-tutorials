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

if(${CMAKE_VERSION} VERSION_GREATER "3.6.3") 
  ctest_configure(CAPTURE_CMAKE_ERROR CONFIGURE_STATUS)
else()
  ctest_configure(RETURN_VALUE CONFIGURE_STATUS)
endif()
message("------> Configure status : ${CONFIGURE_STATUS}")
if(NOT CONFIGURE_STATUS EQUAL 0)
  message(FATAL_ERROR " *** Configure (cmake) process failed *** ")
endif()

# --- Build ---

if(NOT CTEST_BUILD_CONFIGURATION)
  set(CTEST_BUILD_CONFIGURATION "Profiling")
endif()

if(${CMAKE_VERSION} VERSION_GREATER "3.6.3") 
  ctest_build(
      PROJECT_NAME ${current_project}
      CAPTURE_CMAKE_ERROR BUILD_STATUS)
else()
  ctest_build(
      PROJECT_NAME ${current_project}
      RETURN_VALUE BUILD_STATUS)
endif()

message("------> Build status : ${CONFIGURE_STATUS}")
if(NOT BUILD_STATUS EQUAL 0)
  message(FATAL_ERROR " *** Build (cmake) process failed *** ")
endif()


# -- Tests --
if(WITH_TESTS)
  message("---- Start ctest_ctest process ----")
  
  if(${CMAKE_VERSION} VERSION_GREATER "3.6.3") 
    ctest_test(
      PARALLEL_LEVEL NP
      CAPTURE_CMAKE_ERROR TEST_STATUS
      SCHEDULE_RANDOM ON
      RETURN_VALUE TESTS_RESULTS
      )
  else()
    ctest_test(
      PARALLEL_LEVEL NP
      RETURN_VALUE TEST_STATUS
      SCHEDULE_RANDOM ON
      )
  endif()
  message("---- End ctest_ctest process ----")
  message("------> Tests status : ${TEST_STATUS}")
endif()

# -- memory check --
if(CTEST_BUILD_CONFIGURATION MATCHES "Profiling")
  find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)
  set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--quiet --leak-check=full --show-leak-kinds=definite,possible --track-origins=yes --error-limit=no --gen-suppressions=all") 
  set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--quiet --leak-check=full --show-reachable=yes --error-limit=no --gen-suppressions=all") 
  ctest_memcheck(PARALLEL_LEVEL NP QUIET)
endif()

# -- Submission to cdash --
if(DO_SUBMIT)
  message("---- Start ctest_submit process ----")
  ctest_submit(
    CAPTURE_CMAKE_ERROR  SUBMISSION_STATUS
    RETRY_COUNT 4 # Retry 4 times, if submission failed ...)
    RETRY_DELAY 5 # seconds
    )
  message("---- End ctest_submit process ----")
endif()

# ============= Summary =============
message(STATUS "\n============================================ Summary ============================================")
message(STATUS "CTest process for ${current_project} has ended.")
message(STATUS "Ctest model is: ${model}")
message(STATUS "Ctest executed on sources directory : ${CTEST_SOURCE_DIRECTORY}")
message(STATUS "Build name (cdash) : ${CTEST_BUILD_NAME}")
message(STATUS "Site (cdash) : ${CTEST_SITE}")
message(STATUS "=================================================================================================\n")

# tests failed?
if(WITH_TESTS)
  if(NOT TEST_STATUS EQUAL 0 OR NOT TESTS_RESULTS EQUAL 0)
    message(FATAL_ERROR " *** test failure *** ")
  endif()
endif()
# -- Submission failed? --
if(DO_SUBMIT)
  if(NOT SUBMISSION_STATUS EQUAL 0)
    message(WARNING " *** submission failure *** ")
  endif()
endif()


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

