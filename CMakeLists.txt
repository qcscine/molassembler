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

if(NOT CLANG_TIDY_EXE)
  message(STATUS "No clang-tidy found!")
else()
  message(STATUS "Clang-tidy found: ${CLANG_TIDY_EXE}")
  set(DO_CLANG_TIDY "${CLANG_TIDY_EXE}" "-checks=*,-clang-analyzer-alpha.*")
endif()

# set headers
set(HEADERS_all
  CyclicPolygons.h
)

# sources definition
set(SOURCES_all
  CyclicPolygons.cpp
)

# specify targets
add_library(cyclic_polygons STATIC ${SOURCES_all} ${HEADERS_all})

# Compilation options
if( UNIX )
  add_definitions(-fPIC)
endif()

# Switch debug and release build types
if("${CMAKE_BUILD_TYPE}" STREQUAL "Release") 
  message(STATUS "Selected release build")
  # no debug information, disable assert
  target_compile_definitions(cyclic_polygons PRIVATE NDEBUG)
  target_compile_options(cyclic_polygons PUBLIC -O2)
else() # assume Debug
  message(STATUS "Assume debug build")
  target_compile_options(cyclic_polygons PUBLIC -g)
endif()

set_target_properties(
  cyclic_polygons PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON
  COMPILE_FLAGS ${WARNING_FLAGS}
)

include_directories(SYSTEM
  ${EIGEN3_INCLUDE_DIRS}
  "/home/jan-grimo/PhD/projects/libs/"
)

if(CLANG_TIDY_EXE)
  set_target_properties(
    cyclic_polygons PROPERTIES
    CXX_CLANG_TIDY "${DO_CLANG_TIDY}"
  )
endif()

add_executable(cyclicPolygonTests tests.cpp)
target_compile_options(cyclicPolygonTests PUBLIC -O2 -g)
set_target_properties(
  cyclicPolygonTests PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON
  COMPILE_FLAGS ${WARNING_FLAGS}
)
target_include_directories(cyclicPolygonTests PRIVATE ${BOOST_ROOT})
target_link_libraries(cyclicPolygonTests ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} cyclic_polygons)

add_executable(analysis analysis.cpp)
target_compile_options(analysis PUBLIC -O2 -g)
set_target_properties(
  analysis PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON
  COMPILE_FLAGS ${WARNING_FLAGS}
)
target_include_directories(analysis PRIVATE ${BOOST_ROOT})
target_link_libraries(analysis ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} cyclic_polygons)

# add test
add_test(NAME testAll COMMAND analysis)
