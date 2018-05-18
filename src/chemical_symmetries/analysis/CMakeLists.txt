cmake_minimum_required(VERSION 3.0.2)

project(analysis)

find_package(Boost REQUIRED COMPONENTS unit_test_framework program_options system)
find_package(Eigen3 3.1.1 REQUIRED)

add_executable(writeTransitionDotFiles writeTransitionDotFiles.cpp)
target_compile_options(writeTransitionDotFiles PUBLIC -g -O2 -Wall -Wpedantic -Wextra)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  target_compile_options(writeTransitionDotFiles PRIVATE -fconstexpr-steps=1999999999 -fconstexpr-backtrace-limit=0)
endif()
target_link_libraries(writeTransitionDotFiles PRIVATE Boost::boost chemical_symmetries)

add_executable(writeMappings writeMappings.cpp)
target_compile_options(writeMappings PUBLIC -g -O2 -Wall -Wpedantic -Wextra)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  target_compile_options(writeMappings PRIVATE -fconstexpr-steps=1999999999 -fconstexpr-backtrace-limit=0)
endif()
target_link_libraries(writeMappings
  PRIVATE
    Boost::boost
    Boost::program_options
    Boost::system
    chemical_symmetries
)

set_property(TARGET writeMappings PROPERTY FOLDER "Analysis")
set_property(TARGET writeTransitionDotFiles PROPERTY FOLDER "Analysis")
