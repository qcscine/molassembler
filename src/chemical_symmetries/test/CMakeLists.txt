project(SymmetryTestsProject)

find_package(Boost REQUIRED COMPONENTS unit_test_framework program_options)
find_package(Eigen3 3.1.1 REQUIRED)
find_package(temple 0.1.0 REQUIRED)

add_executable(symmetryTests tests.cpp)
target_compile_options(symmetryTests PUBLIC -g -O2)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  target_compile_options(
    symmetryTests
    PRIVATE -fconstexpr-steps=1999999999 -fconstexpr-backtrace-limit=0
  )
endif()

target_include_directories(
  symmetryTests PRIVATE ${EIGEN3_INCLUDE_DIRS}
)
target_link_libraries(symmetryTests
  PUBLIC
    Boost::boost
    Boost::unit_test_framework
    chemical_symmetries
)

# add test
add_test(NAME SymmetryTests COMMAND symmetryTests)

set_property(TARGET symmetryTests PROPERTY FOLDER "Tests")
