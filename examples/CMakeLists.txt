# Set minimum version for cmake
cmake_minimum_required(VERSION 3.5)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#include(WindowsSiconosSetup) Outdated? To be reviewed

set(tests_timeout 35 CACHE INTERNAL "Limit time for tests (in seconds)")

project(siconos_examples CXX)


# -- Try to find a proper siconos installation --

# Get siconos version required to run the examples
include(SiconosRequiredVersion)

# find_package: look for siconosConfig.cmake file
# --> this set siconos_ROOT_DIR to the directory where siconos is installed
# Remark : we may also use the package registry magic to find siconos:
# see https://cmake.org/cmake/help/git-master/manual/cmake-packages.7.html#creating-packages
find_package(siconos ${SICONOS_REQUIRED_VERSION} REQUIRED)
  
# --- Search siconos package path ---
# --> based on path to script siconos
find_file(SICONOS_SCRIPT
  siconos
  PATHS ${siconos_ROOT_DIR}/bin
  )
   
if(NOT SICONOS_SCRIPT)
  message(FATAL_ERROR "Can not find siconos script. Please set your PATH environment variable properly.")
endif()

# Add cmake files from siconos install
list(APPEND CMAKE_MODULE_PATH ${siconos_ROOT_DIR}/share/siconos/cmake)

# Some usefull macros
include(SiconosTools)

# -- Python API --
# 'HAS_PYTHON_WRAPPER' provided by siconosConfig.cmake
# Check if python API has been generated and look for
# some specific packages if required.
if(HAS_PYTHON_WRAPPER)
  message("Siconos Python API is ON")
  set(PYTHON_EXECUTABLE ${SICONOS_PYTHON_EXECUTABLE})
  include(FindPythonModule)
  find_python_module(matplotlib)
  find_python_module(pyhull)
 if(SICONOS_HAS_OCE)
   find_python_module(OCC REQUIRED)
 endif()
endif()

# -- Set the list of examples to be processed --
# It depends on the components installed on the system

# Minimal set of examples (i.e. assume at least numerics component)
set(EXAMPLES_DIRECTORIES
  friction_contact_solvers
  )


# -- Check which components are available --
string(FIND "${installed_targets}" "kernel" kernel_installed )
string(FIND "${installed_targets}" "control" control_installed )
string(FIND "${installed_targets}" "mechanics" mechanics_installed )
# mechanics with bullet? --> SICONOS_HAS_BULLET = ON
# mechanics or mechanisms with OCE? --> SICONOS_HAS_OCE = ON
string(FIND "${installed_targets}" "mechanisms" mechanisms_installed )
string(FIND "${installed_targets}" "io" io_installed )

if(${io_installed} GREATER -1)
  find_python_module(h5py REQUIRED)
endif()


# --- Sort examples, depending on the minimal component they need ---

# Warning : when a directory is listed in one of the variables below
# all its subdirectories are also taken into account!

# It's a bit boring to set explicitely all these dirs. We should review examples organisation ...

# numerics only
set(numerics_examples
  friction_contact_solvers
  )

# kernel
set(kernel_examples
  biology
  electronics
  robotics
  mechanics/BallOnMovingPlane
  mechanics/BouncingBall
  mechanics/BouncingBallEventDriven
  mechanics/CamFollower
  mechanics/ChainsOfBalls
  mechanics/ColumnOfBeads
  mechanics/FEM/ImpactingBar
  mechanics/FEM/exemple_Houari
  mechanics/JointsTestsWithBoundaryConditions
  mechanics/JointsTestsWithInternalForces
  mechanics/Music
  mechanics/NewtonEuler
  mechanics/SliderCrank
  mechanics/TwoBouncingBalls
  mechanics/Woodpecker
  mechanics/yoyo
  )

# Check mfem and add this example
if(HAS_MFEM) # To be done!
  list(APPEND kernel_examples mechanics/FEM/Punch)
endif()

# control
set(control_examples
  control
  )

# mechanics
set(mechanics_examples
  )

# mechanics + bullet
set(mechanics_bullet_examples
  mechanics/BulletBouncingBox
  mechanics/JointsTestsWithContactDetection
  mechanics/JointsTests
  mechanics/Queries
  mechanics/RockFall
  mechanics/InterestingDynamics
  mechanics/RigidBodiesCollection
  mechanics/GeometricsPrimitives
  mechanics/GranularMaterial
  )

# mechanics + oce
set(mechanics_oce_examples
  mechanics/OCC_Examples
  )

# mechanisms # Warning : special care for mechanisms
# Please set explicitely each directory that must be handled.
set(mechanisms_examples
  mechanics/Mechanisms/ShaftTube
  mechanics/Mechanisms/SliderCrank
  mechanics/Mechanisms/Trip_Magnetic
  )

# io

set(io_examples
  )

if(QGLVIEWER_FOUND)
  list(APPEND io_examples mechanics/Spheres/NewtonEuler)
  list(APPEND io_examples mechanics/Spheres/Lagrangian)
endif()


# We assume that numerics, at least, is installed ...
set(EXAMPLES_DIRECTORIES ${numerics_examples})

if(${kernel_installed} GREATER -1)
  message(STATUS "kernel component found")
  set(EXAMPLES_DIRECTORIES  "${EXAMPLES_DIRECTORIES}" ${kernel_examples})
endif()

if(${control_installed} GREATER -1)
  message(STATUS "control component found")
  set(EXAMPLES_DIRECTORIES  "${EXAMPLES_DIRECTORIES}" ${control_examples})
endif()

if(${mechanics_installed} GREATER -1)
  message(STATUS "mechanics component found")
  set(EXAMPLES_DIRECTORIES  "${EXAMPLES_DIRECTORIES}" ${mechanics_examples})
  if(SICONOS_HAS_BULLET)
    message(STATUS "add examples requiring bullet")
    set(EXAMPLES_DIRECTORIES  "${EXAMPLES_DIRECTORIES}" ${mechanics_bullet_examples})
  endif()
  if(SICONOS_HAS_OCE)
    message(STATUS "add examples requiring oce")
    set(EXAMPLES_DIRECTORIES  "${EXAMPLES_DIRECTORIES}" ${mechanics_oce_examples})
  endif()
  if(${io_installed} GREATER -1)
    message(STATUS "io component found")
    set(EXAMPLES_DIRECTORIES  "${EXAMPLES_DIRECTORIES}" ${io_examples})
  endif()
endif()


# Exclude some tests if LpSolve has not been found
if(NOT HAS_ONE_LP_SOLVER)
  set(NO_TEST_FILES ${NO_TEST_FILES} SMCElectroPneumaticItw Twisting)
endif()

# sharing a double between an executable and a library is hackisch
if(MSVC)
  set(NO_TEST_FILES ${NO_TEST_FILES} RockingBlock Pendulum)
endif(MSVC)

#  At this point, two lists are available :
#  - EXAMPLES_DIRECTORIES : all dirs that shoulb be scanned to find examples sources
#  - NO_TEST_FILES : a list of directories that will be excluded from tests.

include(CTest)

if(BUILD_TESTING)
  # -- Scan each dir in examples_directories and collect source files
  foreach(_D ${EXAMPLES_DIRECTORIES})
    message("Scan examples directory ${_D}")
    get_subdirectories(dirlist ${CMAKE_CURRENT_SOURCE_DIR}/${_D})
    if(CMAKE_VERSION VERSION_LESS 3.12.3)
      set(newlist)
      foreach(f ${dirlist})
        list(APPEND newlist "${_D}/${f}")
      endforeach()
      set(dirlist "${newlist}")
    else()
      list(TRANSFORM dirlist PREPEND ${_D}/)
    endif()

    list(APPEND dirlist ${_D})
    if(CMAKE_VERSION VERSION_LESS 3.6)
      list_filter(dirlist src)
      list_filter(dirlist lugin)
    else()
      list(FILTER dirlist EXCLUDE REGEX "src")
      list(FILTER dirlist EXCLUDE REGEX "Plugin")
      list(FILTER dirlist EXCLUDE REGEX "plugin")
    endif()
    foreach(_dir ${dirlist})
      set(src_dir ${CMAKE_CURRENT_SOURCE_DIR}/${_dir})
      if(HAS_PYTHON_WRAPPER)
        file(GLOB EXAMPLES_P ${src_dir}/*.cpp ${src_dir}/*.py)
      else()
        file(GLOB EXAMPLES_P ${src_dir}/*.cpp)
      endif()
      file(GLOB EXAMPLES_DATA ${src_dir}/*.ref ${src_dir}/*/*.mat)
      if(CMAKE_VERSION VERSION_LESS 3.6)
        list_filter(EXAMPLES_P src)
        list_filter(EXAMPLES_P lugin)
      else()
        list(FILTER EXAMPLES_P EXCLUDE REGEX "Plugin")
        list(FILTER EXAMPLES_P EXCLUDE REGEX "plugin")
      endif()
      if(EXAMPLES_P)
        foreach(_P ${EXAMPLES_P})
	  # Full path to current file ...
	  get_filename_component(_RP ${_P} ABSOLUTE)
	  set(TEST_ME TRUE)
	  # Exclude 'no test' files
	  foreach(_NT ${NO_TEST_FILES})
	    if(_RP MATCHES "${_NT}")
	      set(TEST_ME FALSE)
	    endif()
	  endforeach()
	  if(TEST_ME)
	    get_filename_component(EXAMPLE_DIR ${_RP} PATH)
	    get_filename_component(EXAMPLE_NAME ${_RP} NAME_WE)
	    get_filename_component(ext ${_P} EXT)
	    # message("Add ${_D}/${_dir}/${EXAMPLE_NAME}")
            
            set(bin_dir ${CMAKE_CURRENT_BINARY_DIR}/${_dir}/${EXAMPLE_NAME})
	    file(MAKE_DIRECTORY ${bin_dir})
	    foreach(datafile ${EXAMPLES_DATA})
	      #message(STATUS "configure " ${datafile})
	      configure_file(${datafile} ${bin_dir} COPYONLY)
	    endforeach()
            string(REPLACE "/" "_" target_name ${_dir})
            set(target_name ${target_name}_${EXAMPLE_NAME})
            if(ext MATCHES ".py")
 	      set(target_name ${target_name}_py)
 	    else()
	      add_custom_target(${target_name}_build ALL ${CMAKE_COMMAND} -E chdir ${bin_dir} ${SICONOS_SCRIPT} --noexec --silent ${_P}
		WORKING_DIRECTORY ${bin_dir})
	    endif()
            
            # Create a 'test' for the current example
            # name : driver file (cpp or py)  name
            #
 	    add_test(
              NAME ${target_name}
              COMMAND ${CMAKE_COMMAND} -E chdir ${bin_dir} ${SICONOS_SCRIPT} ${_P})
            # set timeout value
	    set_tests_properties(${target_name} PROPERTIES TIMEOUT ${tests_timeout})
          endif()
        endforeach()
      endif()
    endforeach()
  endforeach()

  # separate treatment for mechanisms
  if(${mechanisms_installed} GREATER -1)
    message(STATUS "mechanisms component found")
    foreach(_dir ${mechanisms_examples})
      file(COPY ${_dir} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/mechanisms_examples)
      get_filename_component(target_name ${_dir} NAME)
      set(bin_dir ${CMAKE_CURRENT_BINARY_DIR}/mechanisms_examples/${target_name})
      
      add_test(
        NAME ${target_name}
        COMMAND ${CMAKE_COMMAND} -E chdir ${bin_dir} ${siconos_ROOT_DIR}/bin/siconos_mechanisms .)
      # set timeout value
      set_tests_properties(${target_name} PROPERTIES TIMEOUT ${tests_timeout})
      
    endforeach()
  endif()

  
endif()  

# The following should work, but isn't ... so I used NO_TEST_FILES -- xhub
#SET_TESTS_PROPERTIES("SliderCrankD1MinusLinear" PROPERTIES WILL_FAIL TRUE)


message("============================================================================================")
message("Found siconos script : ${SICONOS_SCRIPT}")
get_filename_component(SICONOS_INSTALL_DIR ${SICONOS_SCRIPT} PATH)
get_filename_component(SICONOS_INSTALL_DIR ${SICONOS_INSTALL_DIR} PATH)
message("The examples will be run with siconos installed in : ${SICONOS_INSTALL_DIR}")
message("Try ctest -N to get a list of all tests that will be processed.")
message("============================================================================================")
