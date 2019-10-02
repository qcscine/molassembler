cmake_minimum_required(VERSION 3.9)

project(molassembler)

# Version information
set(molassembler_VERSION_MAJOR 0)
set(molassembler_VERSION_MINOR 1)
set(molassembler_VERSION_PATCH 0)
set(molassembler_VERSION ${molassembler_VERSION_MAJOR}.${molassembler_VERSION_MINOR}.${molassembler_VERSION_PATCH})
set(molassembler_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# Project-specific options
option(MOLASSEMBLER_ANALYSIS_BINARIES "Compile analysis binaries" OFF)
option(MOLASSEMBLER_VALIDATION_TESTS "Compile validation tests" OFF)
option(MOLASSEMBLER_DOWNLOAD_DEPENDENCIES "Enable downloading dependencies (for CI-Only!)" OFF)
option(MOLASSEMBLER_PARALLELIZE "Enables OpenMP parallelization" ON)
option(MOLASSEMBLER_IPO "Try to enable interprocedural optimization" OFF)
option(MOLASSEMBLER_NO_GNU_UNIQUE "Set --no-gnu-unique GCC flag" OFF)
option(MOLASSEMBLER_ADDRESS_SANITIZE "Add address sanitizer" OFF)

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
