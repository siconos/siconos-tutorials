#[=======================================================================[.rst:
ctest_tools.cmake
----
This file contains some functions used to setup ctest config

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
  cmake_parse_arguments(run "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  message("------> status/result : ${_STATUS}/${_RESULT}")
  if(NOT _STATUS EQUAL 0 OR NOT _RESULT EQUAL 0)
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

  if(${CMAKE_VERSION} VERSION_GREATER "3.10.3") 
    # https://cmake.org/cmake/help/latest/command/cmake_host_system_information.html
    cmake_host_system_information(RESULT osname QUERY OS_NAME)
    cmake_host_system_information(RESULT osplatform QUERY OS_PLATFORM)
    cmake_host_system_information(RESULT hostname QUERY HOSTNAME)
  else()
    set(osname ${CMAKE_SYSTEM_NAME})
    set(osplatform ${CMAKE_SYSTEM_PROCESSOR})
  endif()

  string(STRIP ${osname} osname)
  string(STRIP ${osplatform} osplatform)

  if(CI_GITLAB)
    string(FIND $ENV{CI_JOB_IMAGE} siconos length)
    string(SUBSTRING $ENV{CI_JOB_IMAGE} ${length}+7 -1 dockerimagename)
    string(STRIP ${dockerimagename} dockerimagename)
    set(hostname "[ ${dockerimagename} from gitlab registry]")
  elseif(CI_TRAVIS)
    set(hostname "[ ${hostname} hosted on travis]") 
  endif()
  
  set(_SITE "${osname}-${osplatform}${hostname}")
  string(STRIP _SITE ${_SITE})
  set(CTEST_SITE "${_SITE}" PARENT_SCOPE)
endfunction()

# set build name, according to host, ci, git status ...
function(set_cdash_build_name)
  # Get hash for commit of current version of Siconos
  # Saved in the file siconos-commit.txt during Siconos installation
  # (only if WITH_GIT=ON).
  file(READ ${SICONOS_INSTALL_DIR}/share/siconos-commit.txt branch_commit)
  string(STRIP ${branch_commit} branch_commit)
  include(${CTEST_SOURCE_DIRECTORY}/cmake/SiconosRequiredVersion.cmake)
  set(_name "Examples [based on Siconos-${SICONOS_REQUIRED_VERSION},${branch_commit}]")
  string(STRIP ${_name} _name)
  set(CTEST_BUILD_NAME "${_name}" PARENT_SCOPE)
endfunction()

