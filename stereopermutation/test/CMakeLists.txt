cmake_minimum_required(VERSION 3.0.2)

project(tests)

find_package(Boost REQUIRED COMPONENTS unit_test_framework program_options)
find_package(Eigen3 3.1.1 REQUIRED)

add_executable(permutationTests tests.cpp)
target_compile_options(permutationTests PUBLIC -O2 -g -Wall -Wpedantic -Wextra)
target_include_directories(permutationTests PRIVATE ${BOOST_INCLUDE_DIRS})
target_link_libraries(permutationTests PUBLIC ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} stereopermutation)

# add test target
add_test(NAME test1 COMMAND permutationTests)

SET_PROPERTY(TARGET permutationTests PROPERTY FOLDER "Tests")
