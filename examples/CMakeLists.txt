# Set minimum version for cmake
cmake_minimum_required(VERSION 3.0.2)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Force out-of-source build
include(OutOfSourcesBuild)

# Some usefull macros
include(SiconosTools)

#include(WindowsSiconosSetup)
set(tests_timeout 35 CACHE INTERNAL "Limit time for tests (in seconds)")

project(siconos_examples CXX)

#include(CXXCompilerSetup)
#include(CCompilerSetup)

# -- List of directories to exclude from tests --
set(NO_TEST_FILES QGL tools)

# sharing a double between an executable and a library is hackisch
if(MSVC)
  set(NO_TEST_FILES ${NO_TEST_FILES} RockingBlock Pendulum)
endif(MSVC)


# -- Try to find a proper siconos installation --

# Get siconos current version
include(SiconosVersion)
# find_package: look for siconosConfig.cmake file in CMAKE_MODULE_PATH
# --> this set siconos_ROOT_DIR to the directory where siconos is installed
# Remark : we may also use the package registry magic to find siconos:
# see https://cmake.org/cmake/help/git-master/manual/cmake-packages.7.html#creating-packages
find_package(siconos ${SICONOS_VERSION} REQUIRED)
  
# --- Search siconos package path ---
# --> based on path to script siconos
find_file(SICONOS_SCRIPT
  siconos
  PATHS ${siconos_ROOT_DIR}/bin
  )
   
  
if(NOT SICONOS_SCRIPT)
  message(FATAL_ERROR "Can not find siconos script. Please set your PATH environment variable properly.")
endif()

# -- Set the list of examples to be processed --
# It depends on the components installed on the system

# Minimal set of examples (i.e. assume at least numerics component)
set(EXAMPLES_DIRECTORIES
  friction_contact_solvers
  )

string(FIND "${installed_targets}" "kernel" kernel_installed )
string(FIND "${installed_targets}" "control" control_installed )
string(FIND "${installed_targets}" "io" io_installed )

  
# - Kernel component examples -
if(${kernel_installed} GREATER -1)
  MESSAGE(STATUS "kernel component found")
  set(EXAMPLES_DIRECTORIES
    "${EXAMPLES_DIRECTORIES}"
    biology
    electronics
    mechanics
    robotics/HuMAns_pa10
    #Spheres/Lagrangian
    #Spheres/NewtonEuler
    )
endif()

# - Control component examples -
if(${control_installed} GREATER -1)
  message(STATUS "control component found")
  set(EXAMPLES_DIRECTORIES
    "${EXAMPLES_DIRECTORIES}"
    control
    )
endif()

  
# Some examples need mechanics component. Do not include
# them if mechanics not found.
# Note FP : we should move all examples that require mechanics component in a specific directory
string(FIND "${installed_targets}" "mechanics" mechanics_installed )

# set the list of directories of examples based on mechanics component.
set(mechanics_examples
  JointsTests
  RockFall
  RigidBodiesCollection
  OCC_Examples
  InterestingDynamics
  GranularMaterial
  GeometricPrimitives
  DirectProjection
  Hacks
  Queries
  )

if (${mechanics_installed} LESS 0)
  message("mechanics component not found")
  foreach(_dir ${mechanics_examples})
    list(APPEND NO_TEST_FILES ${_dir})
  endforeach()
endif()

# -- Python API --
# Check if python API has been generated and look for
# some specific packages if required.
if(HAS_PYTHON_WRAPPER)
  set(PYTHON_EXECUTABLE ${SICONOS_PYTHON_EXECUTABLE})
  include(FindPackageHandleStandardArgs)
  include(FindPythonModule)
  find_python_module(matplotlib)
  find_python_module(h5py)
  if(NOT h5py_FOUND)
    list(APPEND NO_TEST_FILES "MechanicsIO")
  endif()
else()
  list(APPEND NO_TEST_FILES "MechanicsIO")
endif()

# Exclude some tests if LpSolve has not been found
if(NOT HAS_ONE_LP_SOLVER)
  set(NO_TEST_FILES ${NO_TEST_FILES} SMCElectroPneumaticItw Twisting)
endif()

# Bullet examples
if(HAS_BULLET)
  list(APPEND EXAMPLES_DIRECTORIES Mechanics/ContactDetection)
  
  if(HAS_PYTHON_WRAPPER)
    if(NOT h5py_FOUND)
      list(APPEND NO_TEST_FILES "BulletIO")
    endif()
  else()
    list(APPEND NO_TEST_FILES "BulletIO")
  endif()
else()
  list(APPEND NO_TEST_FILES BulletBouncingBox)
endif()

#  At this point, two lists are available :
#  - EXAMPLES_DIRECTORIES : all dirs that shoulb be scanned to find examples sources
#  - NO_TEST_FILES : a list of directories that will be excluded from tests.

include(CTest)

if(BUILD_TESTING)
  # -- Scan each dir in examples_directories and collect source files
  foreach(_D ${EXAMPLES_DIRECTORIES})
    message("Scan examples directory ${_D}")
    get_subdirectories(dirlist ${CMAKE_CURRENT_SOURCE_DIR}/${_D})
    list(APPEND dirlist .)
    foreach(_dir ${dirlist})
      set(src_dir ${CMAKE_CURRENT_SOURCE_DIR}/${_D}/${_dir})
      if(WITH_PYTHON_WRAPPER)
        file(GLOB EXAMPLES_P ${src_dir}/*.cpp ${src_dir}/*.py)
      else()
        file(GLOB EXAMPLES_P ${src_dir}/*.cpp)
      endif()
      file(GLOB EXAMPLES_DATA ${src_dir}/*.ref ${src_dir}/*/*.mat)
      
      if(EXAMPLES_P)
        foreach(_P ${EXAMPLES_P})
	  # Full path to current file ...
	  get_filename_component(_RP ${_P} ABSOLUTE)
	  if(NOT _RP MATCHES Plugin) # To exclude plugins
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
	      message("Add ${_D}/${_dir}/${EXAMPLE_NAME}")
              
              set(bin_dir ${CMAKE_CURRENT_BINARY_DIR}/${_D}/${_dir}/${EXAMPLE_NAME})
	      file(MAKE_DIRECTORY ${bin_dir})
	      foreach(datafile ${EXAMPLES_DATA})
	        message(STATUS "configure " ${datafile})
	        configure_file(${datafile} ${bin_dir})
	      endforeach()
	      if(ext MATCHES ".py")
 	        set(target_name ${_dir}_${EXAMPLE_NAME}_py)
 	      else()
	        set(target_name ${_dir}_${EXAMPLE_NAME})
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
	  endif()
        endforeach()
      endif()
    endforeach()
  endforeach()
endif()  

# The following should work, but isn't ... so I used NO_TEST_FILES -- xhub
#SET_TESTS_PROPERTIES("SliderCrankD1MinusLinear" PROPERTIES WILL_FAIL TRUE)


message("============================================================================================")
message("Found siconos script : ${SICONOS_SCRIPT}")
get_filename_component(SICONOS_INSTALL_DIR ${SICONOS_SCRIPT} PATH)
get_filename_component(SICONOS_INSTALL_DIR ${SICONOS_INSTALL_DIR} PATH)
message("The examples will be run with siconos installed in : ${SICONOS_INSTALL_DIR}")
message("============================================================================================")
