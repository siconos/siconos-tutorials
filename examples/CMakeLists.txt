# Set minimum version for cmake
cmake_minimum_required(VERSION 3.19)

cmake_policy(SET CMP0074 NEW) # find_package() uses <PackageName>_ROOT variables.
cmake_policy(SET CMP0086 NEW) # UseSWIG honors SWIG_MODULE_NAME via -module flag.
cmake_policy(SET CMP0079 NEW) # target_link_libraries() allows use with targets in other directories
cmake_policy(SET CMP0067 NEW) # Honor language standard in try_compile() source-file signature.

project(siconos_examples CXX C)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(tests_timeout 300 CACHE INTERNAL "Limit time for tests (in seconds)")

# -- Try to find a proper siconos installation --
# Get siconos version required to run the examples
include(SiconosRequiredVersion)

# find_package looks for a proper installation of siconos and especially for the
# package configuration file siconos-config.cmake.
# See https://cmake.org/cmake/help/git-master/command/find_package.html?highlight=find_package
#
# If siconos is not installed in a standard place, cmake searches in directories listed in CMAKE_PREFIX_PATH.
# This variable must contains the path to siconos installation.
#
# Remark : we may also use the package registry magic to find siconos:
# see https://cmake.org/cmake/help/git-master/manual/cmake-packages.7.html#creating-packages
find_package(siconos REQUIRED)

# --- Search siconos package path ---
# --> based on path to script siconos.
# Note FP : this would probably be better to set this variable in siconos-config.cmake.
# find_file(siconos_SCRIPT
#   siconos
#   NO_DEFAULT_PATH
#   PATHS ${siconos_ROOT}/bin
#   )
if(NOT siconos_SCRIPT)
  message(FATAL_ERROR "Can not find siconos script. Please set your PATH environment variable properly.")
endif()

# Some usefull macros
include(SiconosTools)

# -- Python API --
# 'HAS_PYTHON_WRAPPER' provided by siconosConfig.cmake
# Check if python API has been generated and look for
# some specific packages if required.
if(HAS_PYTHON_WRAPPER OR HAS_PYB11_WRAPPER)
  message("Siconos Python API is ON")
  include(FindPythonModule)
  find_python_module(matplotlib)
  if(SICONOS_HAS_OpenCASCADE)
    find_python_module(OCC REQUIRED)
  endif()
endif()

# --- Sort examples, depending on the minimal component they need ---

# Warning : when a directory is listed in one of the variables below
# all its subdirectories are also taken into account!

# It's a bit boring to set explicitely all these dirs. We should review examples organisation ...

# -- Set the list of examples to be processed --
# It depends on the components installed on the system


#  === numerics ===
# We assume that numerics, at least, is installed ...
set(numerics_examples
  friction_contact_solvers
)
list(APPEND EXAMPLES_DIRECTORIES ${numerics_examples})

# === kernel ===
string(FIND "${installed_targets}" "kernel" kernel_installed)
if(${kernel_installed} GREATER -1)
  message(STATUS "kernel component found")
  set(kernel_examples
    biology/StepSystem
    electronics/BuckConverter
    electronics/CircuitRLCD
    electronics/DiodeBridge
    electronics/DiodeBridgePowSup
    electronics/PowerConverter
    electronics/Circuit-ThetaGammaSchemes
    electronics/ColpittsOscillator
    electronics/DiodeBridgeCapFilter
    electronics/IdealSwitch
    mechanics/2DRigidBodies
    mechanics/BallOnMovingPlane
    mechanics/BouncingBall
    mechanics/BouncingBallEventDriven
    mechanics/CamFollower
    mechanics/ChainsOfBalls
    mechanics/ColumnOfBeads
    mechanics/FEM/ImpactingBar
    mechanics/FEM/native
    mechanics/FourBarLinkageWith3clearance
    mechanics/FrictionOscillator
    mechanics/Masonry
    mechanics/Music
    mechanics/NewtonEuler
    mechanics/Pendulum
    mechanics/PendulumOnASlider
    mechanics/Relay
    mechanics/RockingBlock
    mechanics/SliderCrank
    mechanics/TwoBouncingBalls
    mechanics/Woodpecker
    mechanics/yoyo

  )
  # Check mfem and add this example
  if(HAS_MFEM) # To be done!
    list(APPEND kernel_examples mechanics/FEM/Punch)
  endif()
  list(APPEND EXAMPLES_DIRECTORIES ${kernel_examples})


  # sharing a double between an executable and a library is hackisch
  if(MSVC)
    list(APPEND EXCLUDE_PATTERN_FROM_TESTS RockingBlock Pendulum)
  endif()

endif()

# === vkernel ===
# the vkernel examples sit beside the kernel ones they are compared with
# (mechanics/BouncingBall: *Vkernel*.cpp)
string(FIND "${installed_targets}" "vkernel" vkernel_installed)
if(${vkernel_installed} GREATER -1)
  message(STATUS "vkernel component found")
else()
  list(APPEND EXCLUDE_PATTERN_FROM_TESTS Vkernel)
endif()


# === control ===
string(FIND "${installed_targets}" "control" control_installed)
if(${control_installed} GREATER -1)
  message(STATUS "control component found")
  set(control_examples
    control
  )
  list(APPEND EXAMPLES_DIRECTORIES ${control_examples})
  # Exclude some tests if LpSolve has not been found
  if(NOT HAS_ONE_LP_SOLVER)
    list(APPEND EXCLUDE_PATTERN_FROM_TESTS SMCElectroPneumaticItw Twisting)
  endif()
endif()

# === mechanics ===
string(FIND "${installed_targets}" "mechanics" mechanics_installed)
# mechanics with bullet? --> SICONOS_HAS_BULLET = ON
# mechanics or mechanisms with OpenCASCADE? --> SICONOS_HAS_OpenCASCADE  = ON
# mechanics + bullet
if(${mechanics_installed} GREATER -1)
  message(STATUS "mechanics component found")
  set(mechanics_examples
    mechanics/JointsTestsWithBoundaryConditions
    mechanics/JointsTestsWithInternalForces
  )
  list(APPEND EXAMPLES_DIRECTORIES ${mechanics_examples})
  include(FindPythonModule)

  list(APPEND EXCLUDE_PATTERN_FROM_TESTS meshio_prepost matplotlib_post_2d)

  if(SICONOS_HAS_BULLET)
    message(STATUS "add examples requiring bullet")

    set(mechanics_bullet_examples
      mechanics/BulletBouncingBox
      mechanics/DirectProjection
      mechanics/GeometricPrimitives
      mechanics/GranularMaterial
      mechanics/JointsTestsWithContactDetection
      mechanics/InterestingDynamics
      mechanics/JointsTests
      mechanics/Queries
      mechanics/RockFall
      mechanics/RigidBodiesCollection
      mechanics/RollingFriction
      mechanics/TimeofBirthandDeath
    )
    list(APPEND EXAMPLES_DIRECTORIES ${mechanics_bullet_examples})

    if(HAS_PYB11_WRAPPER)
      find_python_module(pylmgc90)
    endif()
    if(NOT pylmgc90_FOUND)
      list(APPEND EXCLUDE_PATTERN_FROM_TESTS mkspheres)
    endif()

    # RockFall examples depending on optional python packages
    if(HAS_PYTHON_WRAPPER OR HAS_PYB11_WRAPPER)
      find_python_module(stl)      # numpy-stl
      find_python_module(trimesh)
    endif()
    if(NOT stl_file)
      list(APPEND EXCLUDE_PATTERN_FROM_TESTS NSCD_model_adapted_energy_1702_FB)
    endif()
    if(NOT trimesh_file)
      list(APPEND EXCLUDE_PATTERN_FROM_TESTS simu_rockfall_biplan)
    endif()

  endif()
  if(SICONOS_HAS_OpenCASCADE)
    message(STATUS "add examples requiring occt")
    # mechanics + occ
    set(mechanics_occt_examples
      mechanics/MultiBodySystems/OCC_examples
      mechanics/MultiBodySystems/SliderCrank
      mechanics/MultiBodySystems/WatchEscapment
    )
    list(APPEND EXAMPLES_DIRECTORIES ${mechanics_occt_examples})
  endif()
endif()


# === io ===
string(FIND "${installed_targets}" "io" io_installed)
if(${io_installed} GREATER -1)
  message(STATUS "io component found")
  find_python_module(h5py REQUIRED)

  # no examples ?

endif()


#  At this point, two lists are available :
#  - EXAMPLES_DIRECTORIES : all dirs that shoulb be scanned to find examples sources
#  - EXCLUDE_PATTERN_FROM_TESTS : a list of patterns, tests name matching with them will be excluded.
enable_testing()

# -- Scan each dir in examples_directories and collect source files
foreach(current_dir IN LISTS EXAMPLES_DIRECTORIES)
  message("Scan examples directory ${current_dir}")

  # -- collect all source files for the current example dirs --
  # Each file found will be used to create a test.
  set(src_dir ${CMAKE_CURRENT_SOURCE_DIR}/${current_dir})
  if(HAS_PYTHON_WRAPPER OR HAS_PYB11_WRAPPER)
    file(GLOB current_example_files ${src_dir}/*.cpp ${src_dir}/*.py)
  else()
    file(GLOB current_example_files ${src_dir}/*.cpp)
  endif()

  # exclude plugins (they are treated separatly)
  list(FILTER current_example_files EXCLUDE REGEX "Plugin")
  list(FILTER current_example_files EXCLUDE REGEX "plugin")

  # Exclude 'no test' files
  foreach(pattern IN LISTS EXCLUDE_PATTERN_FROM_TESTS)
    list(FILTER current_example_files EXCLUDE REGEX "${pattern}")
  endforeach()

  if(NOT ${mechanics_installed} GREATER -1)
    list(FILTER current_example_files EXCLUDE REGEX "mechanics/2DRigidBodies/Ball2D.cpp")
  endif()

  # Create a working directory for the current example.
  set(working_dir ${CMAKE_CURRENT_BINARY_DIR}/${current_dir})
  file(MAKE_DIRECTORY ${working_dir})

  # collect data files for the current example dirs
  file(GLOB current_example_data ${src_dir}/*.ref ${src_dir}/*/*.mat
    ${src_dir}/*.tess ${src_dir}/coor*.txt ${src_dir}/rad*.txt)
  file(GLOB_RECURSE gmsh_files RELATIVE ${src_dir} ${src_dir}/*.msh*)

  foreach(datafile IN LISTS gmsh_files)
    configure_file(${src_dir}/${datafile} ${working_dir}/${datafile} COPYONLY)
  endforeach()


  # Copy data files into the working directory.
  foreach(datafile IN LISTS current_example_data)
    configure_file(${datafile} ${working_dir} COPYONLY)
  endforeach()

  file(GLOB current_example_dir_data ${src_dir}/data/*.asc)
  # Copy data files into the working directory.
  foreach(datafile IN LISTS current_example_dir_data)
    file(MAKE_DIRECTORY ${working_dir}/data)
    configure_file(${datafile} ${working_dir}/data/ COPYONLY)
  endforeach()

  string(REPLACE "/" "_" test_dir_name ${current_dir})
  # # Create a target to clean all build directories
  # # (.siconos for each _dir)
  # Not needed if we use tmpdir for build
  # add_custom_target(clean_${test_dir_name} ALL
  #   COMMAND ${siconos_SCRIPT} --clean
  #   WORKING_DIRECTORY ${working_dir}
  #   COMMENT "Clean ${working_dir}")

  # -- For each collected source file, create a test --
  # The tests run with SICONOS_EXAMPLES_SHORT_RUN set: the longest examples read it to
  # run a shorter simulation (and compare it with their -short.ref file, if any).
  # Run directly, without this variable, they do the full simulation.
  foreach(driver IN LISTS current_example_files)
    # Set test name (based on example dir name)
    get_filename_component(filename_we ${driver} NAME_WE)
    set(test_name ${test_dir_name}_${filename_we})
    get_filename_component(ext ${driver} EXT)

    # Create test
    if(ext MATCHES ".py")
      # -- Python --
      # Build everything in temp dirs (autom. removed at the end of the test)
      add_test(
        NAME py_${test_name}
        COMMAND ${siconos_SCRIPT} ${driver}
        WORKING_DIRECTORY ${working_dir}
      )

      # python drivers of a directory share its .siconos/ and addons/
      # (built on the fly by the siconos script): never run them concurrently.
      set_tests_properties(py_${test_name} PROPERTIES
        TIMEOUT ${tests_timeout}
        RESOURCE_LOCK ${test_dir_name}
        ENVIRONMENT "SICONOS_EXAMPLES_SHORT_RUN=1")
    else()
      # -- C++ --
      # A test for the build part (siconos --no-exec)
      # Build everything in temp dirs (autom. removed at the end of the test)
      add_test(NAME build_${test_name}
        COMMAND ${siconos_SCRIPT} --no-exec ${driver} --tmp-dir
        WORKING_DIRECTORY ${working_dir})

      # and a test for the run part
      add_test(
        NAME run_${test_name}
        COMMAND ${siconos_SCRIPT} ${driver} --no-build
        WORKING_DIRECTORY ${working_dir})

      set_tests_properties(run_${test_name} PROPERTIES
        DEPENDS build_${test_name}
        TIMEOUT ${tests_timeout}
        ENVIRONMENT "SICONOS_EXAMPLES_SHORT_RUN=1")

    endif()
  endforeach()
endforeach()

# separate treatment for mechanisms
# === mechanisms ===
string(FIND "${installed_targets}" "mechanisms" mechanisms_installed)
if(${mechanisms_installed} GREATER -1)
  message(STATUS "mechanisms component found")

  set(mechanisms_examples
    mechanics/Mechanisms/ShaftTube
    mechanics/Mechanisms/SliderCrank
    mechanics/Mechanisms/Trip_Magnetic
  )
  #list(APPEND EXAMPLES_DIRECTORIES ${mechanisms_examples})
  foreach(current_dir IN LISTS mechanisms_examples)
    set(working_dir ${CMAKE_CURRENT_BINARY_DIR}/${current_dir})
    # file(MAKE_DIRECTORY ${working_dir})

    file(COPY ${current_dir} DESTINATION ${working_dir}) #CMAKE_CURRENT_BINARY_DIR}/mechanisms_examples)
    string(REPLACE "/" "_" test_name ${current_dir})

    add_test(
      NAME ${test_name}
      COMMAND ${siconos_mechanisms_SCRIPT} .
      WORKING_DIRECTORY ${working_dir})
    # set timeout value
    set_tests_properties(${target_name} PROPERTIES TIMEOUT ${tests_timeout})

  endforeach()
endif()

message("============================================================================================")
message("Found siconos : ")
execute_process(COMMAND ${siconos_SCRIPT} --info)

message("\nTry ctest -N to get a list of all tests that will be processed.")
message("============================================================================================")
