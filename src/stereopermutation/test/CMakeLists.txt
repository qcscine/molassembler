project(tests)

if(NOT TARGET Boost::unit_test_framework OR NOT TARGET Boost::program_options)
  find_package(Boost REQUIRED COMPONENTS unit_test_framework program_options)
endif()
find_package(Eigen3 3.1.1 REQUIRED)

add_executable(permutationTests tests.cpp)
target_compile_options(permutationTests PRIVATE ${MOLASSEMBLER_CXX_FLAGS})
target_link_libraries(permutationTests
  PRIVATE
    Boost::boost
    Boost::unit_test_framework
    stereopermutation
)

# add test target
add_test(NAME StereopermutationTests COMMAND permutationTests)

set_property(TARGET permutationTests PROPERTY FOLDER "Tests")
