cmake_minimum_required(VERSION 3.9)

project(molassembler
  VERSION 0.1.0
  DESCRIPTION "Graph-based molecular toolkit for possibly inorganic molecules"
  LANGUAGES CXX
)
enable_testing()

# Project-specific options
option(MOLASSEMBLER_ANALYSIS_BINARIES "Compile analysis binaries" OFF)
option(MOLASSEMBLER_VALIDATION_TESTS "Compile validation tests" OFF)
option(MOLASSEMBLER_PARALLELIZE "Enables OpenMP parallelization" ON)
option(MOLASSEMBLER_IPO "Try to enable interprocedural optimization" OFF)
option(MOLASSEMBLER_NO_GNU_UNIQUE "Set --no-gnu-unique GCC flag" OFF)
option(MOLASSEMBLER_ADDRESS_SANITIZE "Add address sanitizer" OFF)
option(MOLASSEMBLER_BUILD_DOCUMENTATION "Builds and installs documentation" OFF)

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

add_subdirectory(doc)
add_subdirectory(src)
