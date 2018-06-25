project(tests)

set(Boost_USE_STATIC_LIBS OFF)
find_package(Boost REQUIRED COMPONENTS unit_test_framework program_options)
find_package(Eigen3 3.1.1 REQUIRED)

add_executable(permutationTests tests.cpp)
target_compile_options(permutationTests PUBLIC -O2 -g)
target_link_libraries(permutationTests
  PRIVATE
    Boost::boost
    Boost::unit_test_framework
    stereopermutation
)

# add test target
add_test(NAME StereopermutationTests COMMAND permutationTests)

set_property(TARGET permutationTests PROPERTY FOLDER "Tests")
