#-----------------------------------------------
# Ctest driver to build and tests examples.
# Target : continuous integration on gitlab-ci or Travis,
#
# Input variables :
# - SICONOS_INSTALL_DIR : path to Siconos install.
# - OSNAME : host system name (used to qualify cdash build). If not set, try to catch info
#   using common commands (lsb_release ...)
# ----------------------------------------------

# Assumes :
# - a proper install of siconos in SICONOS_INSTALL_DIR
# - an environment variable CI_PROJECT_DIR which contains the path to
# siconos-tutorial repository
# - a file siconos-commit-number.txt in $HOME.
#   to siconos-tutorials repository.


# -- job : build and test examples -- 
message("--- Start conf for siconos examples build and tests.")

# ============= setup  ================

# -- CI_PROJECT_DIR is a required environment variable --
# --> set by default for gitlab-ci, even inside the docker container
if(DEFINED ENV{TRAVIS})
  if($ENV{TRAVIS} STREQUAL true)
    set(CI_TRAVIS ON)
    set(ENV{CI_PROJECT_DIR} ${CTEST_SOURCE_DIRECTORY})
  endif()
endif()
if(DEFINED ENV{GITLAB_CI})
  if($ENV{GITLAB_CI} STREQUAL true)
    set(CI_GITLAB ON)
  endif()
endif()
if(NOT DEFINED ENV{CI_PROJECT_DIR} )
  message(FATAL_ERROR "Please set env variable CI_PROJECT_DIR to siconos-tutorials sources directory (git repo).")
endif()

# -- Definition of all variables required for ctest --
include($ENV{CI_PROJECT_DIR}/ci_gitlab/ctest_tools.cmake)

# ------------------
# Here starts ctest config
# ------------------

# - Source dir and path to siconos install
if(NOT CTEST_SOURCE_DIRECTORY)
  set(CTEST_SOURCE_DIRECTORY $ENV{CI_PROJECT_DIR}/examples)
endif()

# - Top level build directory -
# If not specified : current dir.
if(NOT CTEST_BINARY_DIRECTORY)
  set(CTEST_BINARY_DIRECTORY .)
endif()


# Build name (for cdash)
if(NOT CTEST_BUILD_NAME)
  set_cdash_build_name()
endif()

if(NOT CTEST_CMAKE_GENERATOR)
  set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
endif()

if(NOT CTEST_SITE)
  set_site_name()
endif()

if(NOT CTEST_BUILD_CONFIGURATION)
  set(CTEST_BUILD_CONFIGURATION "Release")
endif()


# =============  Run ctest steps ================

# Current testing model. Priority: 
# Nightly -> set by scheduler on gricad-gitlab
# Continuous -> set in .gitlab-ci.yml
# Experimental : default
if(NOT model)
  set(model Experimental)
endif()

ctest_start(${model})
# Set CTEST_CONFIGURE_COMMAND to cmake followed by siconos options
if(NOT SICONOS_INSTALL_DIR)
  set(SICONOS_INSTALL_DIR /home/install-siconos)
endif()
# setup to find a specific version of siconos.
include($ENV{CI_PROJECT_DIR}/examples/cmake/SiconosRequiredVersion.cmake)
set(ConfigPackageLocation lib/cmake/siconos-${SICONOS_REQUIRED_VERSION})
set(CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND} -Dsiconos_DIR=${SICONOS_INSTALL_DIR}/${ConfigPackageLocation}")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${CTEST_SOURCE_DIRECTORY}")
set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE ${SICONOS_INSTALL_DIR}/${ConfigPackageLocation})

# Remark: do not used parallel for examples. It leads to: warning: jobserver unavailable: using -j1. Add `+' to parent make rule.

#set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE ${CTEST_SOURCE_DIRECTORY}/cmake/valgrind.supp)

message("\n\n=============== Start ctest_configure =============== ")
message("- Configure command line :\n ${CTEST_CONFIGURE_COMMAND}\n")

ctest_configure(
  RETURN_VALUE _RESULT
  CAPTURE_CMAKE_ERROR _STATUS
  #QUIET
  )
post_ctest(PHASE Configure)
message("\n\n=============== Start ctest_build =============== ")

cmake_host_system_information(RESULT NP QUERY NUMBER_OF_LOGICAL_CORES)
if(NOT ALLOW_PARALLEL_BUILD)
  set(NP 1)
endif()
ctest_build(
  FLAGS -j${NP}
  CAPTURE_CMAKE_ERROR _STATUS
  RETURN_VALUE _RESULT
  #QUIET if quiet, travis failed because of missing outputs during a long time ...
  )
post_ctest(PHASE Build)
message("\n\n=============== Start ctest_test (nbprocs = ${NP}) =============== ")
ctest_test(
  #PARALLEL_LEVEL NP
  CAPTURE_CMAKE_ERROR _STATUS
  #SCHEDULE_RANDOM ON
  RETURN_VALUE _RESULT
  # QUIET
  )

if(CDASH_SUBMIT)
  ctest_submit(
    RETURN_VALUE RETURN_STATUS
    CAPTURE_CMAKE_ERROR SUBMISSION_STATUS
    )
  if(NOT SUBMISSION_STATUS EQUAL 0)
    message(WARNING " *** submission failure *** ")
  endif()
endif()
# ============= Summary =============
message(STATUS "\n============================================ Summary ============================================")
message(STATUS "CTest process for siconos-tutorials has ended.")
message(STATUS "Ctest model is: ${model}")
message(STATUS "Ctest executed on sources directory : ${CTEST_SOURCE_DIRECTORY}")
message(STATUS "CTEST_BINARY_DIRECTORY is: ${CTEST_BINARY_DIRECTORY}")
message(STATUS "CTEST_BUILD_CONFIGURATION is: ${CTEST_BUILD_CONFIGURATION}")
if(CDASH_SUBMIT)
  message(STATUS "Cdash server name: ${CTEST_DROP_SITE}/${CTEST_DROP_LOCATION}.")
  message(STATUS "Cdash build name: ${CTEST_BUILD_NAME}")
  message(STATUS "Cdash Site name: ${CTEST_SITE}")
endif()
message(STATUS "=================================================================================================\n")
