cmake_minimum_required(VERSION 3.0.2)

project(tests)

find_package(Boost REQUIRED COMPONENTS unit_test_framework)
find_package(Eigen3 3.1.1 REQUIRED)

# TODO compile_flags is deprecated, introduce compile_options instead

add_executable(templeTests tests.cpp)
target_compile_options(templeTests PUBLIC -g -Og ${WARNING_FLAGS})
target_include_directories(templeTests PRIVATE ${BOOST_ROOT})
target_link_libraries(templeTests PUBLIC temple ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

add_executable(templeConstexprTests constexprTests.cpp)
target_compile_options(templeConstexprTests PUBLIC -g -Og ${WARNING_FLAGS})
target_include_directories(templeConstexprTests PRIVATE ${BOOST_ROOT})
target_link_libraries(templeConstexprTests PUBLIC temple ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# add test
add_test(NAME Base COMMAND templeTests)
add_test(NAME Constexpr COMMAND templeConstexprTests)

SET_PROPERTY(TARGET templeTests PROPERTY FOLDER "Tests")
SET_PROPERTY(TARGET templeConstexprTests PROPERTY FOLDER "Tests")
