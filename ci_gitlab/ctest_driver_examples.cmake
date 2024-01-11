#-----------------------------------------------
# Ctest driver to build and tests examples.
# Target : continuous integration on gitlab-ci or Travis,
#
# Input variables :
# - OSNAME : host system name (used to qualify cdash build). If not set, try to catch info
#   using common commands (lsb_release ...)
# ----------------------------------------------

# Assumes :
# - a proper install of siconos
# - a file siconos-commit-number.txt in $HOME.

# -- job : build and test examples -- 
message("\n\n============== Start conf for siconos examples build and tests ==============\n\n")

# ============= setup  ================

if(NOT CTEST_CMAKE_GENERATOR)
  set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
endif()

# Warning cmake/ctest run on examples dir!
include(${CTEST_SOURCE_DIRECTORY}/../ci_gitlab/ctest_tools.cmake)

# Build name (for cdash)
if(NOT CTEST_BUILD_NAME)
  set_cdash_build_name()
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
include(${CTEST_SOURCE_DIRECTORY}/cmake/SiconosRequiredVersion.cmake)
set(ConfigPackageLocation lib/cmake/siconos-${SICONOS_REQUIRED_VERSION})
set(CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND} ${CTEST_SOURCE_DIRECTORY}")
set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE ${siconos_DIR})

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
#if(NOT ALLOW_PARALLEL_BUILD)
#  set(NP 1)
#endif()
ctest_build(
  FLAGS -j${NP}
  CAPTURE_CMAKE_ERROR _STATUS
  RETURN_VALUE _RESULT
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
post_ctest(PHASE Test FORCE)

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
