# oldest tested version of cmake used (too lazy to figure out newest compatible
# cmake version)
cmake_minimum_required(VERSION 3.5.1)

# little variables
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# declaration
project(GORProject)

set(GOR1_VERSION 0.1.0)

set(HEADERS_All ${CMAKE_CURRENT_SOURCE_DIR}/gor1/Gor1.h)

enable_testing()

add_library(gor1 INTERFACE)
target_sources(
  gor1 INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/gor1/Gor1.h>
)
target_include_directories(
  gor1 INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_include_directories(
  gor1 SYSTEM INTERFACE $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

set(ConfigPackageLocation lib/cmake/gor1)
set(IncludeInstallDir include/)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/gor1Config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/gor1/gor1Config.cmake"
  INSTALL_DESTINATION ${ConfigPackageLocation}
  PATH_VARS IncludeInstallDir
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/gor1/gor1ConfigVersion.cmake"
  VERSION ${GOR1_VERSION}
  COMPATIBILITY AnyNewerVersion
)

install(
  TARGETS gor1
  EXPORT gor1Targets
  DESTINATION include
)

install(
  EXPORT gor1Targets
  DESTINATION "${ConfigPackageLocation}"
)

install(
  FILES ${HEADERS_ALL}
  DESTINATION IncludeInstallDir
)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/gor1/gor1Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/gor1/gor1ConfigVersion.cmake"
  DESTINATION ${ConfigPackageLocation}
)

source_group("Headers" FILES ${HEADERS_All})
