# oldest tested version of cmake used (too lazy to figure out newest compatible
# cmake version)
cmake_minimum_required(VERSION 3.5.1)

# little variables
set(WARNING_FLAGS "-Wall -Wpedantic -Wextra")

# local paths
set(libs_INCLUDE_DIR "/home/jan-grimo/PhD/projects/libs")
set(BOOST_INCLUDE_DIRS "/usr/include/boost")
set(symmetry_LIB_DIR "/home/jan-grimo/PhD/projects/libs/symmetry_information/build")

# declaration
project(StericUniqueness)
enable_testing()

# dependencies
find_package(Boost REQUIRED COMPONENTS unit_test_framework REQUIRED)
find_package(Eigen3 3.1.1 REQUIRED)
find_library(symmetry_LIB symmetries ${symmetry_LIB_DIR})
find_program(
  CLANG_TIDY_EXE
  NAMES "clang-tidy"
  DOC "Path to clang-tidy executable"
)

if(NOT CLANG_TIDY_EXE)
  message(STATUS "No clang-tidy found!")
else()
  message(STATUS "Clang-tidy found: ${CLANG_TIDY_EXE}")
  set(DO_CLANG_TIDY "${CLANG_TIDY_EXE}" "-checks=*,-clang-analyzer-alpha.*")
endif()

include_directories(SYSTEM
  ${EIGEN3_INCLUDE_DIRS}
  ${libs_INCLUDE_DIR}
)

set(HEADERS_all
  Assignment.h
  GenerateUniques.h
  LogicalOperatorTests.h
  Util.h
)

# sources definition
set(SOURCES_all
  Assignment.cpp
  GenerateUniques.cpp
)


# Library target
add_library(Assignment STATIC ${SOURCES_all} ${HEADERS_all})

if( UNIX )
  add_definitions(-fPIC)
endif()

target_compile_options(Assignment PUBLIC -O2 -g)
set_target_properties(
  Assignment PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON
  COMPILE_FLAGS ${WARNING_FLAGS}
)

if(CLANG_TIDY_EXE)
  set_target_properties(
    Assignment PROPERTIES
    CXX_CLANG_TIDY "{DO_CLANG_TIDY}"
  )
endif()

target_link_libraries(Assignment ${symmetry_LIB})


# testing target
add_executable(AssignmentTests tests.cpp)
target_compile_options(AssignmentTests PUBLIC -O2 -g)
set_target_properties(
  AssignmentTests PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON
  COMPILE_FLAGS ${WARNING_FLAGS}
)
if(CLANG_TIDY_EXE)
  set_target_properties(
    AssignmentTests PROPERTIES
    CXX_CLANG_TIDY "{DO_CLANG_TIDY}"
  )
endif()
target_include_directories(AssignmentTests PRIVATE ${BOOST_INCLUDE_DIRS})
target_link_libraries(AssignmentTests ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} Assignment)

# add test target
add_test(NAME test1 COMMAND AssignmentTests)

# properties target
add_executable(Properties properties.cpp)
target_compile_options(Properties PUBLIC -O2 -std=c++14 -Wall -Wpedantic -Wextra)
target_link_libraries(Properties Assignment)
