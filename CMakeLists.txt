# How to make this into a library?
# Dynamic / static? It would be shit if this add-on to Delib would supersede it

# oldest tested version of cmake used (too lazy to figure out newest compatible
# cmake version)
cmake_minimum_required(VERSION 3.5.1)

# local path configuration
set(Delib_LIBRARY_DIR "/home/jan-grimo/PhD/code/delib2-0/build")
set(Delib_INCLUDE_DIR "/home/jan-grimo/PhD/code/delib2-0/src")
set(Assignmentlib_INCLUDE_DIR "/home/jan-grimo/PhD/projects/steric_uniqueness")

# declaration
project(Molecule)

# dependencies
find_package (Eigen3 REQUIRED)
find_library(DELIB_LIBRARY DelibStatic ${Delib_LIBRARY_DIR})

include_directories(SYSTEM
  ${EIGEN3_INCLUDE_DIRS}
  ${Delib_INCLUDE_DIR}
  ${Assignmentlib_INCLUDE_DIR}
  include
)

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
  include/CN4Stereocenter.h
  include/Cache.h
  include/CommonTrig.h
  include/EZStereocenter.h
  include/Edge.h
  include/EdgeList.h
  include/GraphDistanceMatrix.h
  include/IO.h
  include/Molecule.h
  include/StdlibTypeAlgorithms.h
  include/Stereocenter.h
  include/StereocenterList.h
  include/Tree.h
  include/TreeAlgorithms.h
  include/Version.h
  include/common_typedefs.h
  include/numeric.h
)

# sources definition
set(SOURCES_all
  source/DistanceGeometry/DistanceBoundsMatrix.cpp
  source/DistanceGeometry/MetricMatrix.cpp
  source/AdjacencyMatrix.cpp
  source/AtomInfo.cpp
  source/BondDistance.cpp
  source/CN4Stereocenter.cpp
  source/GraphDistanceMatrix.cpp
  source/Molecule.cpp
  source/Stereocenter.cpp
)

# specify targets
add_library(molecule STATIC ${SOURCES_all} ${HEADERS_all})

# Compilation options
if( UNIX )
  add_definitions(-fPIC)
endif()
target_compile_options(molecule PUBLIC -g -std=c++14 -Wall -Wpedantic -Wextra)

# Linking
target_link_libraries(molecule ${DELIB_LIBRARY})
