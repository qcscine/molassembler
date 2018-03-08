# oldest tested version of cmake used (too lazy to figure out newest compatible
# cmake version)
cmake_minimum_required(VERSION 3.5.1)

# little variables
set(WARNING_FLAGS "-Wall -Wpedantic -Wextra")

# declaration
project(SymmetryInformation)
enable_testing()

message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "Compiler version: ${CMAKE_CXX_COMPILER_VERSION}")

# dependencies
find_package(Boost REQUIRED COMPONENTS unit_test_framework)
find_package(Eigen3 3.1.1 REQUIRED)
find_program(
  CLANG_TIDY_EXE
  NAMES "clang-tidy"
  DOC "Path to clang-tidy executable"
)

include_directories(SYSTEM
  ${EIGEN3_INCLUDE_DIRS}
  "/home/jan-grimo/PhD/projects/libs/"
)

add_executable(Tests tests.cpp)
target_compile_options(Tests PUBLIC -g)
set_target_properties(
  Tests PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON
  COMPILE_FLAGS ${WARNING_FLAGS}
)
target_include_directories(Tests PRIVATE ${BOOST_ROOT})
target_link_libraries(Tests ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
