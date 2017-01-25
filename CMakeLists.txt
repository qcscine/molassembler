# oldest tested version of cmake used (too lazy to figure out newest compatible
# cmake version)
cmake_minimum_required(VERSION 3.5.1)

# local paths
set(libs_INCLUDE_DIR "/home/jan-grimo/PhD/projects/libs")
set(BOOST_INCLUDE_DIRS "/usr/include/boost")
set(symmetry_LIB "/home/jan-grimo/PhD/projects/libs/symmetry_information/build/libsymmetries.a")

# declaration
project(StericUniqueness)
enable_testing()

# dependencies
find_package(Boost REQUIRED COMPONENTS unit_test_framework)

include_directories(SYSTEM
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

# specify targets
add_library(Assignment STATIC ${SOURCES_all} ${HEADERS_all})

# Compilation options
if( UNIX )
  add_definitions(-fPIC)
endif()
target_compile_options(Assignment PUBLIC -O2 -g -std=c++14 -Wall -Wpedantic -Wextra)
target_link_libraries(Assignment ${symmetry_LIB})

# testing target
add_executable(AssignmentTests tests.cpp)
target_compile_options(AssignmentTests PUBLIC -O2 -g -std=c++14 -Wall -Wpedantic -Wextra)
target_include_directories(AssignmentTests PRIVATE ${BOOST_INCLUDE_DIRS})
target_link_libraries(AssignmentTests ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} Assignment)

# add test
add_test(NAME test1 COMMAND AssignmentTests)
