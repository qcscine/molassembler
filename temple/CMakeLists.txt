# How to make this into a library?

# oldest tested version of cmake used (too lazy to figure out newest compatible
# cmake version)
cmake_minimum_required(VERSION 3.5.1)

# declaration
project(temple)

# Compilation flags
set(WARNING_FLAGS -Wall -Wpedantic -Wextra)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add testing
enable_testing()

# Some messages
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "Compiler version: ${CMAKE_CXX_COMPILER_VERSION}")

# Subdirectories
add_subdirectory(lib)
add_subdirectory(test)
