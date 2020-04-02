# dependencies
if(NOT TARGET Boost::boost)
  find_package(Boost REQUIRED QUIET)
endif()

set(HEADERS_All
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/constexpr/AngleLookup.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/constexpr/CompileTimeOptions.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/constexpr/Data.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/constexpr/Properties.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/ContinuousMeasures.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/CoordinateSystemTransformation.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/Data.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/Diophantine.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/InertialMoments.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/Partitioner.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/PointGroupElements.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/PointGroups.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/Properties.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/PropertyCaching.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/Shapes.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/TauCriteria.h
)

set(SOURCES_All
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/constexpr/Data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/ContinuousMeasures.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/CoordinateSystemTransformation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/Data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/Diophantine.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/InertialMoments.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/Partitioner.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/PointGroupElements.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/Properties.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/PropertyCaching.cpp
)

add_library(shapes STATIC ${SOURCES_All} ${HEADERS_all})
target_include_directories(shapes
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_include_directories(shapes
  SYSTEM PUBLIC $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)
target_include_directories(shapes SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIRS})
target_link_libraries(shapes PUBLIC Boost::boost temple)
target_link_eigen(shapes PUBLIC)

target_compile_options(shapes PRIVATE ${MOLASSEMBLER_CXX_FLAGS})
if(${SHAPES_TRY_CONSTEXPR})
  target_compile_definitions(shapes PUBLIC SHAPES_TRY_CONSTEXPR)
endif()

install(
  TARGETS shapes
  EXPORT shapesTargets
  DESTINATION lib
)

molassembler_install_component_cmake_files(
  COMPONENT shapes
  EXPORT_NAME shapesTargets
)

molassembler_install_component_headers(COMPONENT shapes DIRECTORY shapes)

add_library(molassembler::shapes ALIAS shapes)

# Source groups for Visual Studio
source_group("Headers" FILES ${HEADERS_All})
source_group("Sources" FILES ${SOURCES_All})