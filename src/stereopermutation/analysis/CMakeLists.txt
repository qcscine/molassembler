cmake_minimum_required(VERSION 3.0.2)

project(analysis)

find_package(Boost REQUIRED COMPONENTS program_options)

add_executable(stereopermutationProperties properties.cpp)
target_compile_options(stereopermutationProperties PUBLIC -O2 -std=c++14)
target_link_libraries(stereopermutationProperties
  PRIVATE Boost::program_options stereopermutation
)

set_property(TARGET stereopermutationProperties PROPERTY FOLDER "Analysis")
