# How to make this into a library?
# Dynamic / static? It would be shit if this add-on to Delib would supercede it

# oldest tested version of cmake used (too lazy to figure out newest compatible
# cmake version)
cmake_minimum_required(VERSION 3.5.1)

# little variables
set(WARNING_FLAGS "-Wall -Wpedantic -Wextra")

# declaration
project(SymmetryInformation)
enable_testing()

# dependencies
find_package(Boost REQUIRED COMPONENTS unit_test_framework)
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
  Symmetries.h
)

# sources definition
set(SOURCES_all
  Symmetries.cpp
)

# specify targets
add_library(symmetries STATIC ${SOURCES_all} ${HEADERS_all})

# Compilation options
if( UNIX )
  add_definitions(-fPIC)
endif()
target_compile_options(symmetries PUBLIC -O2 -g)
set_target_properties(
  symmetries PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON
  COMPILE_FLAGS ${WARNING_FLAGS}
)

if(CLANG_TIDY_EXE)
  set_target_properties(
    symmetries PROPERTIES
    CXX_CLANG_TIDY "{DO_CLANG_TIDY}"
  )
endif()

add_executable(symmetryTests tests.cpp)
target_compile_options(symmetryTests PUBLIC -O2 -g)
set_target_properties(
  symmetryTests PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON
  COMPILE_FLAGS ${WARNING_FLAGS}
)
if(CLANG_TIDY_EXE)
  set_target_properties(
    symmetryTests PROPERTIES
    CXX_CLANG_TIDY "{DO_CLANG_TIDY}"
  )
endif()
target_include_directories(symmetryTests PRIVATE ${BOOST_ROOT})
target_link_libraries(symmetryTests ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} symmetries)

# add test
add_test(NAME testAll COMMAND symmetryTests)
