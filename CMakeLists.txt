# How to make this into a library?
# Dynamic / static? It would be shit if this add-on to Delib would supersede it

# oldest tested version of cmake used (too lazy to figure out newest compatible
# cmake version)
cmake_minimum_required(VERSION 3.5.1)

# local path configuration
set(Delib_LIBRARY_DIR "/home/jan-grimo/PhD/code/delib2-0/build")
set(Delib_INCLUDE_DIR "/home/jan-grimo/PhD/code/delib2-0/include")

# declaration
project(Molecule)

# dependencies
find_package (Eigen3 REQUIRED)
find_library(DELIB_LIBRARY delib ${Delib_LIBRARY_DIR})

include_directories(SYSTEM
  ${EIGEN3_INCLUDE_DIRS}
  ${Delib_INCLUDE_DIR}
  include
)

set(HEADERS_all
  include/AdjacencyList.h
  include/AromaticRing.h
  include/BondDistance.h
  include/CommonTrig.h
  include/Edge.h
  include/EdgeList.h
  include/GraphFeature.h
  include/IO.h
  include/Molecule.h
  include/common_typedefs.h
  include/numeric.h
)

# sources definition
set(SOURCES_all
  source/Molecule.cpp
)

# specify targets
add_library(molecule STATIC ${SOURCES_all} ${HEADERS_all})

# Compilation options
if( UNIX )
  add_definitions(-fPIC)
endif()
target_compile_options(molecule PUBLIC -std=c++14 -Wall -Wpedantic -Wextra)

# Linking
target_link_libraries(molecule ${DELIB_LIBRARY})
