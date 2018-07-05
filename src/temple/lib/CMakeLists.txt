set(HEADERS_All
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/constexpr/Array.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/constexpr/Bitset.h
  ${CMAKE_CURRENT_SOURCE_DIR}/temple/constexpr/Bitmask.h
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

# Make the interface an export-name
install(
  TARGETS temple
  EXPORT templeTargets
  DESTINATION lib
)

molassembler_install_component_cmake_files(
  COMPONENT temple
  EXPORT_NAME templeTargets
)

molassembler_install_component_headers(COMPONENT temple DIRECTORY temple)

add_library(molassembler::temple ALIAS temple)

# Source groups for IDEs
source_group("Headers" FILES ${HEADERS_All})
