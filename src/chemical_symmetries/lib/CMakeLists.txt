# dependencies
find_package(Eigen3 3.1.1 REQUIRED)
find_package(Boost REQUIRED)

set(chemical_symmetries_VERSION 0.1.0)

set(HEADERS_All
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/AngleLookup.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/CompileTimeOptions.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/ConstexprProperties.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/DynamicProperties.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Names.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Primitives.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Properties.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Symmetries.h
)

set(SOURCES_All
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/DynamicProperties.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Primitives.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Properties.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Symmetries.cpp
)

add_library(chemical_symmetries STATIC ${SOURCES_All} ${HEADERS_all})
target_include_directories(
  chemical_symmetries
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_include_directories(
  chemical_symmetries
    SYSTEM PUBLIC $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)
target_include_directories(chemical_symmetries PUBLIC ${EIGEN3_INCLUDE_DIRS})
target_link_libraries(chemical_symmetries
  PUBLIC
    Boost::boost
    molassembler::temple
)

target_compile_options(chemical_symmetries PUBLIC)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  target_compile_options(
    chemical_symmetries
    PRIVATE -fconstexpr-steps=1999999999 -fconstexpr-backtrace-limit=0
  )
endif()
if(${CHEMICAL_SYMMETRIES_TRY_CONSTEXPR})
  target_compile_definitions(
    chemical_symmetries PUBLIC CHEMICAL_SYMMETRIES_TRY_CONSTEXPR
  )
endif()

install(
  TARGETS chemical_symmetries
  EXPORT chemical_symmetriesTargets
  DESTINATION lib
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
