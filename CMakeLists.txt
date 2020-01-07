cmake_minimum_required(VERSION 3.9)

project(molassembler
  VERSION 1.0.0
  DESCRIPTION "Graph-based molecular toolkit for possibly inorganic molecules"
  LANGUAGES CXX
)
enable_testing()

# Build molassembler shared by default
option(BUILD_SHARED_LIBS ON)

# Project-specific options
option(MOLASSEMBLER_ANALYSIS_BINARIES "Compile analysis binaries" OFF)
option(MOLASSEMBLER_VALIDATION_TESTS "Compile validation tests" OFF)
option(MOLASSEMBLER_PARALLELIZE "Enables OpenMP parallelization" ON)
option(MOLASSEMBLER_IPO "Try to enable interprocedural optimization" OFF)
option(MOLASSEMBLER_NO_GNU_UNIQUE "Set --no-gnu-unique GCC flag" OFF)
option(MOLASSEMBLER_ADDRESS_SANITIZE "Add address sanitizer" OFF)

# SCINE options
option(SCINE_BUILD_TESTS "Build test executables" ON)
option(SCINE_BUILD_DOCS "Build the documentation" ON)

# Set a default build type
set(default_build_type "Release")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  set(default_build_type "Debug")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to default '${default_build_type}'")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
    STRING "Choose the type of the build."
    FORCE
  )
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY
    STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo"
  )
endif()

if(SCINE_BUILD_DOCS)
  add_subdirectory(doc)
endif()

add_subdirectory(src)
