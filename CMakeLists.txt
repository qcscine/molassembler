# oldest tested version of cmake used (too lazy to figure out newest compatible
# cmake version)
cmake_minimum_required(VERSION 3.5.1)

# declaration
project(chemical_symmetries)

# compilation flags
set(WARNING_FLAGS "-Wall -Wpedantic -Wextra")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Clang constexpr options
set(CLANG_CONSTEXPR_OPTIONS "-fconstexpr-steps=1999999999 -fconstexpr-backtrace-limit=0")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  SET(CONSTEXPR_OPTIONS "${CLANG_CONSTEXPR_OPTIONS}")
endif()

enable_testing()

message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "Compiler version: ${CMAKE_CXX_COMPILER_VERSION}")

add_subdirectory(lib)
add_subdirectory(test)
add_subdirectory(analysis)
