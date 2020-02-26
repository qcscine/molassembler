set(HEADERS_ALL ${CMAKE_CURRENT_SOURCE_DIR}/CyclicPolygons.h)

add_library(cyclic_polygons INTERFACE)
target_sources(cyclic_polygons INTERFACE $<BUILD_INTERFACE:${HEADERS_ALL}>)

target_include_directories(cyclic_polygons INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

target_link_libraries(cyclic_polygons INTERFACE temple)

# dependencies
if(NOT TARGET Boost::unit_test_framework)
  find_package(Boost REQUIRED COMPONENTS unit_test_framework)
endif()

if(SCINE_BUILD_TESTS)
  add_executable(cyclicPolygonTests tests.cpp)

  include(LibType)
  determine_lib_type(Boost::unit_test_framework)
  if(${LIB_TYPE} STREQUAL "SHARED_LIBRARY")
    target_compile_definitions(cyclicPolygonTests PRIVATE BOOST_ALL_DYN_LINK)
  endif()
  unset(LIB_TYPE)

  target_compile_options(cyclicPolygonTests PRIVATE ${MOLASSEMBLER_CXX_FLAGS})
  target_link_libraries(cyclicPolygonTests
    PRIVATE
      Boost::boost
      Boost::unit_test_framework
      temple
  )
  target_link_eigen(cyclicPolygonTests PRIVATE)

  add_test(NAME MASMCyclicPolygonTests COMMAND cyclicPolygonTests)
endif()

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
