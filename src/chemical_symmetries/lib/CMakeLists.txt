# dependencies
find_package(Eigen3 3.1.1 REQUIRED)
if(NOT TARGET Boost::boost)
  find_package(Boost REQUIRED)
endif()

set(chemical_symmetries_VERSION 0.1.0)

set(HEADERS_All
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/AngleLookup.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/CompileTimeOptions.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/ConstexprProperties.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/ContinuousMeasures.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/CoordinateSystemTransformation.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Diophantine.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/DynamicProperties.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/InertialMoments.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Partitioner.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/PointGroupElements.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/PointGroups.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Primitives.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Properties.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Shapes.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Symmetries.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/TauCriteria.h
)

set(SOURCES_All
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/CoordinateSystemTransformation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/ContinuousMeasures.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Diophantine.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/DynamicProperties.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/InertialMoments.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Partitioner.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/PointGroupElements.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Primitives.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Properties.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Symmetries.cpp
)

add_library(chemical_symmetries STATIC ${SOURCES_All} ${HEADERS_all})
target_include_directories(chemical_symmetries
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_include_directories(chemical_symmetries
  SYSTEM PUBLIC $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)
target_include_directories(chemical_symmetries SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIRS})
target_link_libraries(chemical_symmetries
  PUBLIC
    Boost::boost
    temple
)

target_compile_options(chemical_symmetries PRIVATE ${MOLASSEMBLER_CXX_FLAGS})
if(${CHEMICAL_SYMMETRIES_TRY_CONSTEXPR})
  target_compile_definitions(chemical_symmetries
    PUBLIC CHEMICAL_SYMMETRIES_TRY_CONSTEXPR
  )
endif()

# NOTE: The entirety of the generator expression below is provided purely for
# the convenience of installing multiple build types into a single prefix folder
install(
  TARGETS chemical_symmetries
  EXPORT chemical_symmetriesTargets
  DESTINATION $<$<NOT:$<CONFIG:Release>>:$<CONFIG>/>lib
)

molassembler_install_component_cmake_files(
  COMPONENT chemical_symmetries
  EXPORT_NAME chemical_symmetriesTargets
)

molassembler_install_component_headers(
  COMPONENT chemical_symmetries
  DIRECTORY chemical_symmetries
)

add_library(molassembler::chemical_symmetries ALIAS chemical_symmetries)

# Source groups for Visual Studio
source_group("Headers" FILES ${HEADERS_All})
source_group("Sources" FILES ${SOURCES_All})
