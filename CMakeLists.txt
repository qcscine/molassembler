# oldest tested version of cmake used, feel free to lower this if older versions
# work too
cmake_minimum_required(VERSION 3.5.1)

# base paths
set(BASE_DIR "/home/jan-grimo/PhD")
set(LIBS_DIR "${BASE_DIR}/projects/libs")

# local path configuration
set(AssignmentLib_LIBRARY_DIR "${LIBS_DIR}/geometry_assignment/build")
set(Symmetry_LIBRARY_DIR "${LIBS_DIR}/symmetry_information/build")
set(RDL_LIBRARY_DIR "${BASE_DIR}/lib/cpp/RingDecomposerLib/build/src/RingDecomposerLib")
set(RDL_INCLUDE_DIR "${BASE_DIR}/lib/cpp/RingDecomposerLib/src")

# project declaration
project(Molecule)
enable_testing()

message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "Compiler version: ${CMAKE_CXX_COMPILER_VERSION}")

# dependencies: Eigen3
find_package(Eigen3 3.1.1 REQUIRED)
if(Eigen3_FOUND)
  message(STATUS "Eigen version: ${Eigen3_VERSION}")
  message(STATUS "Eigen dir: ${EIGEN3_INCLUDE_DIRS}")
  include_directories(EIGEN3_INCLUDE_DIRS)
endif()

# dependencies: dlib
find_package(dlib 19.4 REQUIRED)
if(dlib_FOUND)
  message(STATUS "dlib version: ${dlib_VERSION}")
  message(STATUS "dlib dir: ${dlib_INCLUDE_DIRS}")
  include_directories(dlib_INCLUDE_DIRS)
endif()

# dependencies: blas (required for dlib)
find_package(BLAS REQUIRED)

# dependencies: Boost
find_package(Boost REQUIRED COMPONENTS unit_test_framework program_options)
if(Boost_FOUND) 
  include_directories(Boost_INCLUDE_DIRS)
endif()

#dependencies: Delib
find_package(Delib REQUIRED)
if(Delib_FOUND)
  message(STATUS "Delib version: ${Delib_VERSION}")
  # unnecessary to add Delib_INCLUDE_DIRS, this is imported through the interface
endif()

# dependencies: libraries
find_library(ASSIGNMENT_LIBRARY Assignment ${AssignmentLib_LIBRARY_DIR})
find_library(SYMMETRY_LIBRARY symmetries ${Symmetry_LIBRARY_DIR})
find_library(RDL_LIBRARY RingDecomposerLibStatic ${RDL_LIBRARY_DIR})

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
  # Enable all clang-tidy checks, deactivating alpha checks and a specific one
  # that trips up over casserts
  set(DO_CLANG_TIDY "${CLANG_TIDY_EXE}" "-checks=*,-clang-analyzer-alpha.*,-cppcoreguidelines-pro-bounds-array-to-pointer-decay")
endif()

# compiler warning flags 
set(WARNING_FLAGS "-Wall -Wpedantic -Wextra")

# import assignment and symmetry libraries so we can specify a dependency
add_library(assignment STATIC IMPORTED)
set_property(TARGET assignment PROPERTY IMPORTED_LOCATION ${ASSIGNMENT_LIBRARY})
add_library(symmetry STATIC IMPORTED)
set_property(TARGET symmetry PROPERTY IMPORTED_LOCATION ${SYMMETRY_LIBRARY})
add_library(rdl STATIC IMPORTED)
set_property(TARGET rdl PROPERTY IMPORTED_LOCATION ${RDL_LIBRARY})
add_dependencies(assignment symmetry)

# add include directories
include_directories(SYSTEM
  ${LIBS_DIR}
  ${RDL_INCLUDE_DIR}
  include
)

# collect headers
set(HEADERS_all
  include/DistanceGeometry/RefinementProblem.h
  include/DistanceGeometry/DistanceBoundsMatrix.h
  include/DistanceGeometry/DistanceGeometry.h
  include/DistanceGeometry/MetricMatrix.h
  include/DistanceGeometry/MoleculeSpatialModel.h
  include/DistanceGeometry/generateConformation.h
  include/AdjacencyMatrix.h
  include/AnalysisHelpers.h
  include/AtomInfo.h
  include/BondDistance.h
  include/CNStereocenter.h
  include/Cache.h
  include/CommonTrig.h
  include/CycleData.h
  include/DelibHelpers.h
  include/EZStereocenter.h
  include/Edges.h
  include/GraphAlgorithms.h
  include/GraphDistanceMatrix.h
  include/IO.h
  include/IterateStereocenterPermutations.h
  include/LocalGeometryModel.h
  include/Molecule.h
  include/MolGraphWriter.h
  include/PersistentIndexMoleculeInterface.h
  include/RangeForTemporary.h
  include/StdlibTypeAlgorithms.h
  include/Stereocenter.h
  include/StereocenterList.h
  include/VSEPR.h
  include/Version.h
  include/common_typedefs.h
)

# collect sources
set(SOURCES_all
  source/DistanceGeometry/DistanceBoundsMatrix.cpp
  source/DistanceGeometry/generateConformation.cpp
  source/DistanceGeometry/MetricMatrix.cpp
  source/DistanceGeometry/MoleculeSpatialModel.cpp
  source/DistanceGeometry/RefinementProblem.cpp
  source/AdjacencyMatrix.cpp
  source/AtomInfo.cpp
  source/BondDistance.cpp
  source/CNStereocenter.cpp
  source/CommonTrig.cpp
  source/CycleData.cpp
  source/DelibHelpers.cpp
  source/EZStereocenter.cpp
  source/GraphDistanceMatrix.cpp
  source/GraphAlgorithms.cpp
  source/LocalGeometryModel.cpp
  source/Log.cpp
  source/Molecule.cpp
  source/MolGraphWriter.cpp
  source/RankingTree.cpp
  source/StdlibTypeAlgorithms.cpp
  source/Stereocenter.cpp
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
  target_compile_options(molecule PUBLIC -O3 -march=native)
else() # assume Debug
  message(STATUS "Assume debug build")
  target_compile_options(molecule PUBLIC -Og -g)
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
    CXX_CLANG_TIDY "${DO_CLANG_TIDY}"
  )
endif()

# Linking
target_link_libraries(molecule Delib::DelibStatic assignment symmetry blas rdl)

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
  if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    target_compile_options(${targetName} PUBLIC -O3 -march=native)
  else()
    target_compile_options(${targetName} PUBLIC -Og -g)
  endif()
  set_target_properties(
    ${targetName} PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON
    COMPILE_FLAGS ${WARNING_FLAGS}
  )

  target_include_directories(${targetName} PRIVATE ${Boost_INCLUDE_DIRS})
  target_link_libraries(
    ${targetName} 
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    Delib::DelibStatic
    assignment
    symmetry
    molecule
    blas
    rdl
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
  if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    target_compile_options(${targetName} PUBLIC -O3 -march=native)
  else()
    target_compile_options(${targetName} PUBLIC -Og -g)
  endif()
  set_target_properties(
    ${targetName} PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON
    COMPILE_FLAGS ${WARNING_FLAGS}
  )

  target_include_directories(${targetName} PRIVATE ${Boost_INCLUDE_DIRS})
  target_link_libraries(
    ${targetName} 
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    Delib::DelibStatic
    assignment
    symmetry
    molecule
    blas
    rdl
  )
endforeach()
