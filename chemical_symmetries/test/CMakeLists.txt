cmake_minimum_required(VERSION 3.0.2)

project(tests)

find_package(Boost REQUIRED COMPONENTS unit_test_framework program_options)
find_package(Eigen3 3.1.1 REQUIRED)
find_package(temple 0.1.0 REQUIRED)

add_executable(symmetryTests tests.cpp)
target_compile_options(symmetryTests PUBLIC -g -O2 -Wall -Wpedantic -Wextra)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  target_compile_options(symmetryTests PRIVATE -fconstexpr-steps=1999999999 -fconstexpr-backtrace-limit=0)
endif()

target_include_directories(symmetryTests PRIVATE ${Boost_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS})
target_link_libraries(symmetryTests 
  PUBLIC 
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} 
    chemicalsymmetries 
)

# add test
add_test(NAME testAll COMMAND symmetryTests)

SET_PROPERTY(TARGET symmetryTests PROPERTY FOLDER "Tests")
