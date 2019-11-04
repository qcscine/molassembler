if(NOT TARGET Boost::unit_test_framework OR NOT TARGET Boost::program_options)
  find_package(Boost REQUIRED COMPONENTS unit_test_framework program_options)
endif()
find_package(Eigen3 3.1.1 REQUIRED)
find_package(temple 0.1.0 REQUIRED)

set(TestFiles
  ${CMAKE_CURRENT_SOURCE_DIR}/ContinuousMeasures.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CoordinateSystemTransformation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Diophantine.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/InertialMoments.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Partitioner.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/tests.cpp
)

add_executable(shapeTests ${TestFiles})

target_include_directories(shapeTests PRIVATE ${EIGEN3_INCLUDE_DIRS})
target_link_libraries(shapeTests
  PUBLIC
    Boost::boost
    Boost::unit_test_framework
    shapes
)

target_compile_options(shapeTests PRIVATE ${MOLASSEMBLER_CXX_FLAGS})
target_compile_definitions(shapeTests PRIVATE BOOST_ALL_DYN_LINK)

# add test
add_test(NAME MASMShapeTests COMMAND shapeTests)

set_property(TARGET shapeTests PROPERTY FOLDER "Tests")
