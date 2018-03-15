cmake_minimum_required(VERSION 3.0.2)

project(analysis)

find_package(Boost REQUIRED COMPONENTS unit_test_framework program_options)

add_executable(stereopermutationProperties properties.cpp)
target_compile_options(stereopermutationProperties PUBLIC -O2 -std=c++14 -Wall -Wpedantic -Wextra)
target_link_libraries(
  stereopermutationProperties
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  stereopermutation
)

SET_PROPERTY(TARGET stereopermutationProperties PROPERTY FOLDER "Analysis")
