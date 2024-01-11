#[=======================================================================[.rst:
ctest_tools.cmake
----
This file contains the functions used to setup ctest config

#]=======================================================================]

#[=======================================================================[.rst:
Call this function after each ctest step (ctest_configure,ctest_build ...)
to handle errors and submission to cdash

Usage :

post_ctest(PHASE <phase_name>)

with phase_name in (Configure, Build, Test).
#]=======================================================================]
function(post_ctest)
  set(oneValueArgs PHASE)
  set(options FORCE)
  cmake_parse_arguments(run "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  message("------> status/result : ${_STATUS}/${_RESULT}")
  # This function is supposed to submit to cdash
  # only if cmake terminates with an error (e.g. FATAL_ERROR) or ctest fails
  if(NOT _STATUS EQUAL 0 OR NOT _RESULT EQUAL 0 OR run_FORCE)
    if(CDASH_SUBMIT)
      ctest_submit(
        RETURN_VALUE RETURN_STATUS
        CAPTURE_CMAKE_ERROR SUBMISSION_STATUS
        )
      if(NOT SUBMISSION_STATUS EQUAL 0)
        message(WARNING " *** submission failure *** ")
      endif()
    else()
      message("- Results won't be submitted to a cdash server.\n")
      return()
    endif()
  endif()
  
  if(NOT _STATUS EQUAL 0 OR NOT _RESULT EQUAL 0)
    message(FATAL_ERROR "\n\n *** ${run_PHASE} process failed *** \n\n")
  endif()
  unset(_RESULT PARENT_SCOPE)
  unset(_STATUS PARENT_SCOPE)
  message("=============== End of ctest ${run_PHASE} =============== ")

endfunction()

# Set site name (cdash) according to current host status. 
function(set_site_name)
  
  # -- Query host system information --
  # --> to set ctest site for cdash.
  #include(cmake_host_system_information)
  cmake_host_system_information(RESULT hostname QUERY HOSTNAME)
  cmake_host_system_information(RESULT fqdn QUERY FQDN)

  cmake_host_system_information(RESULT osname QUERY OS_NAME)
  cmake_host_system_information(RESULT osplatform QUERY OS_PLATFORM)
  cmake_host_system_information(RESULT hostname QUERY HOSTNAME)

  string(STRIP ${osname} osname)
  string(STRIP ${osplatform} osplatform)

  if(DEFINED ENV{GITLAB_CI})
    if($ENV{GITLAB_CI} STREQUAL true)
      string(SUBSTRING $ENV{CI_JOB_IMAGE} 39 -1 dockerimagename) 
      string(STRIP ${dockerimagename} dockerimagename)
      set(hostname "[ ${dockerimagename} from gitlab registry]")
    endif()
  endif()

  set(_SITE "${osname}-${osplatform}-host=${hostname}")

  string(STRIP _SITE ${_SITE})
  set(CTEST_SITE "${_SITE}" PARENT_SCOPE)
endfunction()

# set build name, according to host, ci, git status ...
function(set_cdash_build_name)
  # Get hash for commit of current version of Siconos
  # Saved in the file siconos-commit.txt during Siconos installation
  # (only if WITH_GIT=ON).
  
  execute_process(COMMAND
    siconos --configpath
    OUTPUT_VARIABLE siconos_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  find_file(SICONOS_COMMIT siconos-commit.txt ${siconos_DIR})
  if(SICONOS_COMMIT)
    file(READ ${SICONOS_COMMIT} branch_commit)
    string(STRIP ${branch_commit} branch_commit)
    set(branch_commit "-${branch_commit}")
  endif()
  include(${CTEST_SOURCE_DIRECTORY}/cmake/SiconosRequiredVersion.cmake)
  set(_name "Examples [based on Siconos-${SICONOS_REQUIRED_VERSION}${branch_commit}]")
  string(STRIP ${_name} _name)
  set(CTEST_BUILD_NAME "${_name}" PARENT_SCOPE)
endfunction()

