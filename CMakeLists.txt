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
option(MOLASSEMBLER_DOWNLOAD_DEPENDENICES "Enable downloading dependencies (for CI-Only!)" OFF)
option(MOLASSEMBLER_IPO "Try to enable interprocedural optimization" OFF)
option(MOLASSEMBLER_NO_GNU_UNIQUE "Set --no-gnu-unique GCC flag" OFF)

add_subdirectory(doc)
add_subdirectory(src)
