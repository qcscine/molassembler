# dependencies
find_package(Eigen3 3.1.1 REQUIRED)
if(NOT TARGET Boost::boost)
  find_package(Boost REQUIRED)
endif()

set(HEADERS_All
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/AngleLookup.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/CompileTimeOptions.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/ConstexprProperties.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/ContinuousMeasures.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/CoordinateSystemTransformation.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/Diophantine.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/DynamicProperties.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/InertialMoments.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/Partitioner.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/PointGroupElements.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/PointGroups.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/Primitives.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/Properties.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/Shapes.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/Symmetries.h
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/TauCriteria.h
)

set(SOURCES_All
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/CoordinateSystemTransformation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/ContinuousMeasures.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/Diophantine.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/DynamicProperties.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/InertialMoments.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/Partitioner.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/PointGroupElements.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/Primitives.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/Properties.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/shapes/Symmetries.cpp
)

add_library(shapes STATIC ${SOURCES_All} ${HEADERS_all})
target_include_directories(shapes
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_include_directories(shapes
  SYSTEM PUBLIC $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)
target_include_directories(shapes SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIRS})
target_link_libraries(shapes
  PUBLIC
    Boost::boost
    temple
)

target_compile_options(shapes PRIVATE ${MOLASSEMBLER_CXX_FLAGS})
if(${SHAPES_TRY_CONSTEXPR})
  target_compile_definitions(shapes
    PUBLIC SHAPES_TRY_CONSTEXPR
  )
endif()

# NOTE: The entirety of the generator expression below is provided purely for
# the convenience of installing multiple build types into a single prefix folder
install(
  TARGETS shapes
  EXPORT shapesTargets
  DESTINATION $<$<NOT:$<CONFIG:Release>>:$<CONFIG>/>lib
)

molassembler_install_component_cmake_files(
  COMPONENT shapes
  EXPORT_NAME shapesTargets
)

molassembler_install_component_headers(
  COMPONENT shapes
  DIRECTORY shapes
)

add_library(molassembler::shapes ALIAS shapes)

# Source groups for Visual Studio
source_group("Headers" FILES ${HEADERS_All})
source_group("Sources" FILES ${SOURCES_All})
