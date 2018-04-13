# oldest tested version of cmake used (too lazy to figure out newest compatible
# cmake version)
cmake_minimum_required(VERSION 3.5.1)

# little variables
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# declaration
project(CyclicPolygons)

set(cyclic_polygons_VERSION 0.1.0)

enable_testing()

set(HEADERS_ALL ${CMAKE_CURRENT_SOURCE_DIR}/CyclicPolygons.h)

add_library(cyclic_polygons INTERFACE)
target_sources(
  cyclic_polygons INTERFACE $<BUILD_INTERFACE:${HEADERS_ALL}>
)

target_include_directories(
  cyclic_polygons INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

# dependencies
set(Boost_USE_STATIC_LIBS OFF)
find_package(Boost REQUIRED COMPONENTS unit_test_framework)
find_package(Eigen3 3.1.1 REQUIRED)

add_executable(cyclicPolygonTests tests.cpp)
target_compile_options(cyclicPolygonTests PRIVATE -g -Wall -Wpedantic -Wextra)
target_include_directories(cyclicPolygonTests PRIVATE ${BOOST_ROOT})
target_link_libraries(
  cyclicPolygonTests ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} temple
)

add_test(NAME CyclicPolygonTests COMMAND cyclicPolygonTests)

install(
  TARGETS cyclic_polygons
  EXPORT cyclic_polygonsTargets
  DESTINATION ${INCLUDE_INSTALL_DIR}
)

molassembler_install_component_cmake_files(
  COMPONENT cyclic_polygons
  EXPORT_NAME cyclic_polygonsTargets
)

install(
  FILES ${HEADERS_ALL}
  DESTINATION ${INCLUDE_INSTALL_DIR}
)
