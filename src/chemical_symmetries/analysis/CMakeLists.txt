cmake_minimum_required(VERSION 3.0.2)

project(analysis)

find_package(Boost REQUIRED COMPONENTS unit_test_framework program_options)
find_package(Eigen3 3.1.1 REQUIRED)
find_package(temple 0.1.0 REQUIRED)

add_executable(writeTransitionDotFiles writeTransitionDotFiles.cpp)
target_compile_options(writeTransitionDotFiles PUBLIC -g -O2 -Wall -Wpedantic -Wextra)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  target_compile_options(writeTransitionDotFiles PRIVATE -fconstexpr-steps=1999999999 -fconstexpr-backtrace-limit=0)
endif()
target_include_directories(writeTransitionDotFiles PRIVATE ${Boost_INCLUDE_DIRS} chemicalsymmetries)
target_link_libraries(writeTransitionDotFiles chemicalsymmetries)


message(STATUS "Program options library: ${Boost_PROGRAM_OPTIONS_LIBRARY}")
add_executable(writeMappings writeMappings.cpp)
target_compile_options(writeMappings PUBLIC -g -O2 -Wall -Wpedantic -Wextra)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  target_compile_options(writeMappings PRIVATE -fconstexpr-steps=1999999999 -fconstexpr-backtrace-limit=0)
endif()
target_include_directories(writeMappings PRIVATE ${Boost_INCLUDE_DIRS} chemicalsymmetries)
target_link_libraries(
  writeMappings 
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  chemicalsymmetries
)

SET_PROPERTY(TARGET writeMappings PROPERTY FOLDER "Analysis")
SET_PROPERTY(TARGET writeTransitionDotFiles PROPERTY FOLDER "Analysis")
