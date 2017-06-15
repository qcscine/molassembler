# How to make this into a library?
# Dynamic / static? It would be shit if this add-on to Delib would supercede it

# oldest tested version of cmake used (too lazy to figure out newest compatible
# cmake version)
cmake_minimum_required(VERSION 3.5.1)

# little variables
set(WARNING_FLAGS "-Wall -Wpedantic -Wextra")

message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "Compiler version: ${CMAKE_CXX_COMPILER_VERSION}")

# declaration
project(TemplateMagic)
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


# testing binary
add_executable(magicTests tests.cpp)
target_compile_options(magicTests PUBLIC -O2 -g)
set_target_properties(
  magicTests PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON
  COMPILE_FLAGS ${WARNING_FLAGS}
)
if(CLANG_TIDY_EXE)
  set_target_properties(
    magicTests PROPERTIES
    CXX_CLANG_TIDY "${DO_CLANG_TIDY}"
  )
endif()
target_include_directories(magicTests PRIVATE ${BOOST_ROOT})
target_link_libraries(magicTests ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# add test
add_test(NAME testAll COMMAND magicTests)
