# oldest tested version of cmake used (too lazy to figure out newest compatible
# cmake version)
cmake_minimum_required(VERSION 3.5.1)

# little variables
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# declaration
project(GORProject)

set(gor1_VERSION 0.1.0)

# dependencies
set(Boost_USE_STATIC_LIBS OFF)
find_package(Boost REQUIRED COMPONENTS unit_test_framework)

# Files
set(HEADERS_All ${CMAKE_CURRENT_SOURCE_DIR}/gor1/Gor1.h)

enable_testing()

# Interface library
add_library(gor1 INTERFACE)
target_sources(
  gor1 INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/gor1/Gor1.h>
)
target_include_directories(
  gor1 INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_include_directories(
  gor1 SYSTEM INTERFACE $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

# Installation
install(
  TARGETS gor1
  EXPORT gor1Targets
  DESTINATION ${INCLUDE_INSTALL_DIR}
)

molassembler_install_component_cmake_files(
  COMPONENT gor1
  EXPORT_NAME gor1Targets
)

molassembler_install_component_headers(
  COMPONENT gor1
  DIRECTORY gor1
)

# Tests
add_executable(gor1Tests test.cpp)
if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  target_compile_options(gor1Tests PRIVATE -O3 -march=native)
else()
  target_compile_options(gor1Tests PRIVATE -g)
endif()
target_include_directories(gor1Tests PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(gor1Tests PRIVATE ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

add_test(NAME gor1Tests COMMAND gor1Tests)

set_property(TARGET gor1Tests PROPERTY FOLDER "Tests")

# Source groups
source_group("Headers" FILES ${HEADERS_All})
