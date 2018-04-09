# oldest tested version of cmake used (too lazy to figure out newest compatible
# cmake version)
cmake_minimum_required(VERSION 3.5.1)

# little variables
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# declaration
project(CyclicPolygons)

set(CYCLIC_POLYGONS_VERSION 0.1.0)

enable_testing()

add_library(cyclic_polygons INTERFACE)
target_sources(
  cyclic_polygons INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CyclicPolygons.h>
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

set(ConfigPackageLocation lib/cmake/cyclic_polygons)
set(IncludeInstallDir include/)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cyclic_polygonsConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cyclic_polygons/cyclic_polygonsConfig.cmake"
  INSTALL_DESTINATION ${ConfigPackageLocation}
  PATH_VARS IncludeInstallDir
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/cyclic_polygons/cyclic_polygonsConfigVersion.cmake"
  VERSION ${CYCLIC_POLYGONS_VERSION}
  COMPATIBILITY AnyNewerVersion
)

install(
  TARGETS cyclic_polygons
  EXPORT cyclic_polygonsTargets
  DESTINATION include
)

install(
  EXPORT cyclic_polygonsTargets
  DESTINATION "${ConfigPackageLocation}"
)

install(
  FILES ${HEADERS_ALL}
  DESTINATION IncludeInstallDir
)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/cyclic_polygons/cyclic_polygonsConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/cyclic_polygons/cyclic_polygonsConfigVersion.cmake"
  DESTINATION ${ConfigPackageLocation}
)

source_group("Headers" FILES ${HEADERS_All})
