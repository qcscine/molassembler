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
option(MOLASSEMBLER_ANALYSIS "Compile analysis binaries" OFF)
option(MOLASSEMBLER_VALIDATION "Compile validation tests" OFF)
option(MOLASSEMBLER_PARALLELIZE "Enables OpenMP parallelization" ON)
option(MOLASSEMBLER_IPO "Try to enable interprocedural optimization" OFF)
option(MOLASSEMBLER_NO_GNU_UNIQUE "Set --no-gnu-unique GCC flag" OFF)
option(MOLASSEMBLER_SANITIZE "Add address and UB sanitizers" OFF)

# SCINE options
option(SCINE_BUILD_TESTS "Build test executables" ON)
option(SCINE_BUILD_DOCS "Build the documentation" ON)
option(SCINE_BUILD_PYTHON_BINDINGS "Build python bindings" OFF)

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

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/dev/cmake)
include(ColorMessages)

# Handle sanitization
if(MOLASSEMBLER_SANITIZE)
  cmessage(STATUS "MOLASSEMBLER_SANITIZE: Adding address and UB sanitizers")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address,undefined")
  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -fsanitize=address,undefined")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address,undefined")

  cmessage(WARNING
    "Address sanitiziation can cause ABI incompatibility issues "
    "if dependencies have not been compiled the same way."
  )
endif()

# Handle coverage flag
if(COVERAGE AND CMAKE_BUILD_TYPE STREQUAL "Debug")
  cmessage(STATUS "COVERAGE: Adding coverage flags")
  set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -coverage ${CMAKE_CXX_FLAGS_DEBUG}")
endif()

# Fetch some dependencies already
include(ImportUtilsOS)
import_utils_os()
include(CompilerSpecificSetup)
include(AddEigen)
set(Boost_USE_MULTITHREADED ON)
find_package(Boost
  COMPONENTS regex unit_test_framework program_options filesystem system
  QUIET
)

# Subdirectories
find_package(RingDecomposerLib QUIET)
find_package(nauty QUIET)
add_subdirectory(extern)
add_subdirectory(src)

if(SCINE_BUILD_DOCS)
  include(DoxygenDocumentation)
  include(DoxygenSettings)
  set_doxygen_vars()
  scine_component_documentation(UtilsOSDocumentation)
  unset_doxygen_vars()
endif()
if(SCINE_BUILD_TESTS)
  add_subdirectory(test)
endif()
if(MOLASSEMBLER_ANALYSIS)
  add_subdirectory(analysis)
endif()
if(MOLASSEMBLER_VALIDATION)
  add_subdirectory(validation)
endif()
if(SCINE_BUILD_PYTHON_BINDINGS)
  add_subdirectory(python)
endif()
