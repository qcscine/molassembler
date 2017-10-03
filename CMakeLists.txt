# How to make this into a library?

# oldest tested version of cmake used (too lazy to figure out newest compatible
# cmake version)
cmake_minimum_required(VERSION 3.5.1)

# little variables
set(WARNING_FLAGS "-Wall -Wpedantic -Wextra")
set(CLANG_CONSTEXPR_OPTIONS "-fconstexpr-steps=1999999999 -fconstexpr-backtrace-limit=0")

# declaration
project(SymmetryInformation)
enable_testing()

message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "Compiler version: ${CMAKE_CXX_COMPILER_VERSION}")

# dependencies
find_package(Boost REQUIRED COMPONENTS unit_test_framework)
find_package(Eigen3 3.1.1 REQUIRED)
if(Eigen3_FOUND)
  message(STATUS "Eigen version: ${Eigen3_VERSION}")
  message(STATUS "Eigen dir: ${EIGEN3_INCLUDE_DIRS}")
  include_directories(EIGEN3_INCLUDE_DIRS)
endif()

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
  CompileTimeOptions.h
  ConstexprProperties.h
  DynamicProperties.h
  Primitives.h
  Properties.h
  Symmetries.h
)

# sources definition
set(SOURCES_all
  DynamicProperties.cpp
  Primitives.cpp
  Properties.cpp
  Symmetries.cpp
)

# specify targets
add_library(symmetries STATIC ${SOURCES_all} ${HEADERS_all})

# Compilation options
if( UNIX )
  add_definitions(-fPIC)
endif()

# CLang constexpr options
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  SET(CONSTEXPR_OPTIONS "${CLANG_CONSTEXPR_OPTIONS}")
endif()

# Switch debug and release build types
if("${CMAKE_BUILD_TYPE}" STREQUAL "Release") 
  message(STATUS "Selected release build")
  # no debug information, disable assert
  #  target_compile_definitions(symmetries PRIVATE NDEBUG)
  #target_compile_options(symmetries PUBLIC -O2)
else() # assume Debug
  message(STATUS "Assume debug build")
  target_compile_options(symmetries PUBLIC -g -O2)
endif()

set_target_properties(
  symmetries PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON
  COMPILE_FLAGS "${WARNING_FLAGS} ${CONSTEXPR_OPTIONS}"
)

include_directories(SYSTEM
  ${EIGEN3_INCLUDE_DIRS}
  "/home/jan-grimo/PhD/projects/libs/"
)

#if(CLANG_TIDY_EXE)
#  set_target_properties(
#    symmetries PROPERTIES
#    CXX_CLANG_TIDY "${DO_CLANG_TIDY}"
#  )
#endif()

add_executable(symmetryTests tests.cpp)
target_compile_options(symmetryTests PUBLIC -g -O2)
set_target_properties(
  symmetryTests PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON
  COMPILE_FLAGS "${WARNING_FLAGS} ${CONSTEXPR_OPTIONS}"
)
target_include_directories(symmetryTests PRIVATE ${BOOST_ROOT})
target_link_libraries(symmetryTests ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} symmetries)

# add test
add_test(NAME testAll COMMAND symmetryTests)

add_executable(writeTransitionDotFiles writeTransitionDotFiles.cpp)
target_compile_options(writeTransitionDotFiles PUBLIC -g -O2)
set_target_properties(
  writeTransitionDotFiles PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON
  COMPILE_FLAGS "${WARNING_FLAGS} ${CONSTEXPR_OPTIONS}"
)
target_include_directories(writeTransitionDotFiles PRIVATE ${BOOST_ROOT} symmetries)
target_link_libraries(writeTransitionDotFiles ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} symmetries)
