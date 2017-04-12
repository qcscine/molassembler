# oldest tested version of cmake used, feel free to lower this if older versions
# work too
cmake_minimum_required(VERSION 3.5.1)

# base paths
set(BASE_DIR "/home/jan-grimo/PhD")
set(LIBS_DIR "${BASE_DIR}/projects/libs")

# local path configuration
set(AssignmentLib_LIBRARY_DIR "${LIBS_DIR}/steric_uniqueness/build")
set(Symmetry_LIBRARY_DIR "${LIBS_DIR}/symmetry_information/build")
set(BOOST_INCLUDE_DIRS "/usr/include/boost")
set(DELIB_INCLUDE_DIR "/usr/local/include/Delib")

# project declaration
project(Molecule)
enable_testing()

# dependencies: packages
find_package(Eigen3 3.1.1 REQUIRED)
find_package(Boost REQUIRED COMPONENTS unit_test_framework)
find_package(Delib REQUIRED)
message(STATUS "Delib version: ${Delib_VERSION}")

# dependencies: libraries
find_library(ASSIGNMENT_LIBRARY Assignment ${AssignmentLib_LIBRARY_DIR})
find_library(SYMMETRY_LIBRARY symmetries ${Symmetry_LIBRARY_DIR})

# dependencies: programs
find_program(
  CLANG_TIDY_EXE
  NAMES "clang-tidy"
  DOC "Path to clang-tidy executable"
)
if(NOT CLANG_TIDY_EXE)
  message(STATUS "No clang-tidy found!")
else()
  message(STATUS "Clang-tidy found: ${CLANG_TIDY_EXE}")
  set(DO_CLANG_TIDY "${CLANG_TIDY_EXE}" "-checks=*,-clang-analyzer-alpha.*")
endif()

# compiler warning flags 
set(WARNING_FLAGS "-Wall -Wpedantic -Wextra")

# import assignment and symmetry libraries so we can specify a dependency
add_library(assignment STATIC IMPORTED)
set_property(TARGET assignment PROPERTY IMPORTED_LOCATION ${ASSIGNMENT_LIBRARY})
add_library(symmetry STATIC IMPORTED)
set_property(TARGET symmetry PROPERTY IMPORTED_LOCATION ${SYMMETRY_LIBRARY})
add_dependencies(assignment symmetry)

# add include directories
include_directories(SYSTEM
  ${EIGEN3_INCLUDE_DIRS}
  ${DELIB_INCLUDE_DIR}
  ${LIBS_DIR}
  include
)

# collect headers
set(HEADERS_all
  include/DistanceGeometry/DGRefinementProblem.h
  include/DistanceGeometry/DistanceBoundsMatrix.h
  include/DistanceGeometry/DistanceGeometry.h
  include/DistanceGeometry/MetricMatrix.h
  include/DistanceGeometry/generateConformation.h
  include/AdjacencyList.h
  include/AdjacencyListAlgorithms.h
  include/AdjacencyMatrix.h
  include/AtomInfo.h
  include/BondDistance.h
  include/CNStereocenter.h
  include/Cache.h
  include/CommonTrig.h
  include/EZStereocenter.h
  include/Edges.h
  include/GraphDistanceMatrix.h
  include/IO.h
  include/LocalGeometryModel.h
  include/Molecule.h
  include/PersistentIndexMoleculeInterface.h
  include/StdlibTypeAlgorithms.h
  include/Stereocenter.h
  include/StereocenterList.h
  include/Traits.h
  include/Tree.h
  include/TreeAlgorithms.h
  include/Version.h
  include/VSEPR.h
  include/common_typedefs.h
  include/numeric.h
)

# collect sources
set(SOURCES_all
  source/DistanceGeometry/DistanceBoundsMatrix.cpp
  source/DistanceGeometry/generateConformation.cpp
  source/DistanceGeometry/MetricMatrix.cpp
  source/AdjacencyList.cpp
  source/AdjacencyListAlgorithms.cpp
  source/AdjacencyMatrix.cpp
  source/AtomInfo.cpp
  source/BondDistance.cpp
  source/CNStereocenter.cpp
  source/CommonTrig.cpp
  source/EZStereocenter.cpp
  source/GraphDistanceMatrix.cpp
  source/LocalGeometryModel.cpp
  source/Log.cpp
  source/Molecule.cpp
  source/StdlibTypeAlgorithms.cpp
  source/Stereocenter.cpp
  source/SymmetryFit.cpp
  source/VSEPR.cpp
)

# specify library target
add_library(molecule STATIC ${SOURCES_all} ${HEADERS_all})

# Compilation options
if( UNIX )
  add_definitions(-fPIC)
endif()

# Switch debug and release build types
if("${CMAKE_BUILD_TYPE}" STREQUAL "Release") 
  message(STATUS "Selected release build")
  # no debug information, disable assert
  target_compile_definitions(molecule PRIVATE NDEBUG)
else() # assume Debug
  message(STATUS "Assume debug build")
  target_compile_options(molecule PUBLIC -g)
endif()

# library properties
set_target_properties(
  molecule PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON
  COMPILE_FLAGS ${WARNING_FLAGS}
)

# enable clang-tidy warnings
if(CLANG_TIDY_EXE)
  set_target_properties(
    molecule PROPERTIES
    CXX_CLANG_TIDY "{DO_CLANG_TIDY}"
  )
endif()

# Linking
target_link_libraries(molecule Delib::DelibStatic assignment symmetry)

# Testing targets
# get all test cpp files
file(GLOB testFiles "tests/*.cpp")
foreach(filename ${testFiles})
  # get a target name from the filename
  set(targetName "")
  get_filename_component(targetName ${filename} NAME_WE)
  set(targetName "test_${targetName}")

  # add the executable and specify compilation and linking options
  add_executable(${targetName} ${filename})
  target_compile_options(${targetName} PUBLIC -g)
  set_target_properties(
    ${targetName} PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON
    COMPILE_FLAGS ${WARNING_FLAGS}
  )

  if(CLANG_TIDY_EXE)
    set_target_properties(
      ${targetName} PROPERTIES
      CXX_CLANG_TIDY "{DO_CLANG_TIDY}"
    )
  endif()

  target_include_directories(${targetName} PRIVATE ${BOOST_INCLUDE_DIRS})
  target_link_libraries(
    ${targetName} 
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    Delib::DelibStatic
    assignment
    symmetry
    molecule
  )

  # add target to tests
  add_test(NAME ${targetName} COMMAND ${targetName})
endforeach()

# Analysis targets
file(GLOB analysisFiles "analysis/*.cpp")
foreach(filename ${analysisFiles})
  set(targetName "")
  get_filename_component(targetName ${filename} NAME_WE)
  set(targetName "analysis_${targetName}")

  # add the executable and specify compilation and linking options
  add_executable(${targetName} ${filename})
  target_compile_options(${targetName} PUBLIC -g)
  set_target_properties(
    ${targetName} PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON
    COMPILE_FLAGS ${WARNING_FLAGS}
  )

  if(CLANG_TIDY_EXE)
    set_target_properties(
      ${targetName} PROPERTIES
      CXX_CLANG_TIDY "{DO_CLANG_TIDY}"
    )
  endif()

  target_include_directories(${targetName} PRIVATE ${BOOST_INCLUDE_DIRS})
  target_link_libraries(
    ${targetName} 
    Delib::DelibStatic
    assignment
    symmetry
    molecule
  )
endforeach()
