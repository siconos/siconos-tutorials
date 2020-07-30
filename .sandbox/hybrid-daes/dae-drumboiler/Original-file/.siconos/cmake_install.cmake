# Install script for directory: /scratch/arocca/non-smooth-dae/Jupyter-Python/Siconos-examples/DrumBoiler/DrumBoiler-07-2020/Original-file/.siconos

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/scratch/arocca/non-smooth-dae/Jupyter-Python/Siconos-examples/DrumBoiler/DrumBoiler-07-2020/Original-file/drumboiler-doublemixture" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/scratch/arocca/non-smooth-dae/Jupyter-Python/Siconos-examples/DrumBoiler/DrumBoiler-07-2020/Original-file/drumboiler-doublemixture")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/scratch/arocca/non-smooth-dae/Jupyter-Python/Siconos-examples/DrumBoiler/DrumBoiler-07-2020/Original-file/drumboiler-doublemixture"
         RPATH "/scratch/arocca/.local/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/scratch/arocca/non-smooth-dae/Jupyter-Python/Siconos-examples/DrumBoiler/DrumBoiler-07-2020/Original-file/drumboiler-doublemixture")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/scratch/arocca/non-smooth-dae/Jupyter-Python/Siconos-examples/DrumBoiler/DrumBoiler-07-2020/Original-file" TYPE EXECUTABLE FILES "/scratch/arocca/non-smooth-dae/Jupyter-Python/Siconos-examples/DrumBoiler/DrumBoiler-07-2020/Original-file/.siconos/drumboiler-doublemixture")
  if(EXISTS "$ENV{DESTDIR}/scratch/arocca/non-smooth-dae/Jupyter-Python/Siconos-examples/DrumBoiler/DrumBoiler-07-2020/Original-file/drumboiler-doublemixture" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/scratch/arocca/non-smooth-dae/Jupyter-Python/Siconos-examples/DrumBoiler/DrumBoiler-07-2020/Original-file/drumboiler-doublemixture")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/scratch/arocca/non-smooth-dae/Jupyter-Python/Siconos-examples/DrumBoiler/DrumBoiler-07-2020/Original-file/drumboiler-doublemixture"
         OLD_RPATH "/scratch/arocca/.local/lib:"
         NEW_RPATH "/scratch/arocca/.local/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/scratch/arocca/non-smooth-dae/Jupyter-Python/Siconos-examples/DrumBoiler/DrumBoiler-07-2020/Original-file/drumboiler-doublemixture")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/scratch/arocca/non-smooth-dae/Jupyter-Python/Siconos-examples/DrumBoiler/DrumBoiler-07-2020/Original-file/.siconos/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
