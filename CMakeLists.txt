# How to make this into a library?
# Dynamic / static? It would be shit if this add-on to Delib would supercede it

# oldest tested version of cmake used (too lazy to figure out newest compatible
# cmake version)
cmake_minimum_required(VERSION 3.5.1)

# declaration
project(SymmetryInformation)
enable_testing()

# dependencies
find_package(Boost REQUIRED COMPONENTS unit_test_framework)

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
target_compile_options(symmetries PUBLIC -O2 -g -std=c++14 -Wall -Wpedantic -Wextra)

add_executable(symmetryTests tests.cpp)
target_compile_options(symmetryTests PUBLIC -O2 -g -std=c++14 -Wall -Wpedantic -Wextra)
target_include_directories(symmetryTests PRIVATE ${BOOST_ROOT})
target_link_libraries(symmetryTests ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} symmetries)

# add test
add_test(NAME testAll COMMAND symmetryTests)
