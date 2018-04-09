set(temple_VERSION 0.1.0)

set(HEADERS_All
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/constexpr/Array.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/constexpr/BTree.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/constexpr/ConsecutiveCompare.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/constexpr/Containers.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/constexpr/DynamicArray.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/constexpr/DynamicMap.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/constexpr/DynamicSet.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/constexpr/DynamicUIntArray.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/constexpr/FloatingPointComparison.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/constexpr/LogicalOperatorTests.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/constexpr/Math.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/constexpr/Numeric.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/constexpr/Optional.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/constexpr/Pair.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/constexpr/ToSTL.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/constexpr/TupleType.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/constexpr/TypeTests.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/constexpr/UIntArray.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/constexpr/UpperTriangularMatrix.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/constexpr/Vector.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/Adaptors.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/AddToContainer.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/BTree.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/Cache.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/ContainerTraits.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/Containers.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/Enumerate.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/Functional.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/Invoke.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/MemberFetcher.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/Optionals.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/Random.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/Stringify.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/TinySet.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/Traits.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/UnorderedSets.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/VectorView.h
)

add_library(temple INTERFACE)
target_sources(temple INTERFACE $<BUILD_INTERFACE:${HEADERS_All}>)
target_include_directories(
  temple INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_include_directories(
  temple SYSTEM INTERFACE $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

# Install all relevant header files in <installdir>/include
include(CMakePackageConfigHelpers)

# Set destination for the generated CMake config files
set(ConfigPackageLocation lib/cmake/temple)

set(INCLUDE_INSTALL_DIR include/)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/temple/templeConfig.cmake"
  INSTALL_DESTINATION ${ConfigPackageLocation}
  PATH_VARS INCLUDE_INSTALL_DIR
)

# Generate a basic configuration file in the binary director and subdir temple
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/temple/templeConfigVersion.cmake"
  VERSION ${temple_VERSION}
  COMPATIBILITY AnyNewerVersion
)

install(
  TARGETS temple
  EXPORT templeTargets
  DESTINATION lib
)

install(
  EXPORT templeTargets
  DESTINATION "${ConfigPackageLocation}"
)

install(DIRECTORY
  "${CMAKE_CURRENT_SOURCE_DIR}/temple"
  DESTINATION
  ${INCLUDE_INSTALL_DIR}
  FILES_MATCHING
  PATTERN "*.h"
)

# Install templeConfig.cmake and templeConfigVersion into specified directory
install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/temple/templeConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/temple/templeConfigVersion.cmake"
  DESTINATION ${ConfigPackageLocation}
)

# Source groups for Visual Studio
source_group("Headers" FILES ${HEADERS_All})
