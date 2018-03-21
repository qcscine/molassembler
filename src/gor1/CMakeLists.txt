# oldest tested version of cmake used (too lazy to figure out newest compatible
# cmake version)
cmake_minimum_required(VERSION 3.5.1)

# little variables
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# declaration
project(GORProject)

enable_testing()

add_library(gor1 INTERFACE)
target_sources(gor1 INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/gor1/Gor1.h)
target_include_directories(
  gor1 INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_include_directories(
  gor1 SYSTEM INTERFACE $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)
