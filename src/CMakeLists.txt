# Dependencies
if(MOLASSEMBLER_PARALLELIZE)
  find_package(OpenMP REQUIRED QUIET)
endif()

if(NOT TARGET Boost::filesystem OR NOT TARGET Boost::system)
  find_package(Boost REQUIRED COMPONENTS filesystem system QUIET)
endif()

# Collect headers and sources
file(GLOB_RECURSE HEADERS_all ${CMAKE_CURRENT_SOURCE_DIR}/Molassembler/*.h)
file(GLOB_RECURSE SOURCES_all ${CMAKE_CURRENT_SOURCE_DIR}/Molassembler/*.cpp)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Molassembler/Version.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/Molassembler/Version.h
  @ONLY
)

# Object library to avoid compiling twice for shared and static variants
add_library(molassembler_obj OBJECT ${SOURCES_all})

include(GenerateExportHeader)
generate_export_header(molassembler_obj
  BASE_NAME masm
  EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/Molassembler/Export.h
)

set_target_properties(molassembler_obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(molassembler_obj PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

# Object library build requirements
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.12.0")
  # this is clearly the preferred method of handling this mess, but only
  # available from 3.12 onward. The proper linking an interface population is
  # handled below
  target_link_libraries(molassembler_obj PRIVATE
    Boost::boost
    Scine::UtilsOS
    outcome
    RingDecomposerLibStatic
    json
    Boost::filesystem
    Boost::system
    nauty
    $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
  )
else()
  target_include_directories(molassembler_obj
    PUBLIC
      $<TARGET_PROPERTY:Scine::UtilsOS,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:outcome,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:RingDecomposerLibStatic,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:Boost::boost,INTERFACE_INCLUDE_DIRECTORIES>
  )
endif()
add_eigen(molassembler_obj PRIVATE)
target_compile_options(molassembler_obj PRIVATE
  ${MOLASSEMBLER_CXX_FLAGS}
  $<$<BOOL:${OpenMP_CXX_FOUND}>:${OpenMP_CXX_FLAGS}>
)
# Suppress warnings coming from various external libraries
target_include_directories(molassembler_obj SYSTEM PRIVATE
  ${Boost_INCLUDE_DIR}
  ${RDL_INCLUDE_DIR}
  $<TARGET_PROPERTY:outcome,INTERFACE_INCLUDE_DIRECTORIES>
  $<TARGET_PROPERTY:Eigen3::Eigen,INTERFACE_INCLUDE_DIRECTORIES>
  $<TARGET_PROPERTY:json,INTERFACE_INCLUDE_DIRECTORIES>
  $<TARGET_PROPERTY:nauty,INTERFACE_INCLUDE_DIRECTORIES>
)
set_target_properties(molassembler_obj PROPERTIES
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON
)

function(molassembler_includes_and_properties target_name)
  target_include_directories(${target_name} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
  )
  set_target_properties(${target_name} PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
    OUTPUT_NAME "molassembler"
  )
endfunction()

function(molassembler_library_links target_name)
  # Some of the dependencies of this possibly static version of molassembler
  # are INTERFACE libraries. We need them only for implementation details and
  # they do not provide any symbols, so they are linked PRIVATE. CMake still
  # references their targets as LINK_ONLY dependencies in
  # INTERFACE_LINK_LIBRARIES for the static library!
  #
  # Yet we want only the include information from these libraries so that
  # the library can be built properly, yet not the link dependency.
  # See in particular: https://gitlab.kitware.com/cmake/cmake/issues/15415
  target_include_directories(${target_name} PRIVATE
    $<TARGET_PROPERTY:json,INTERFACE_INCLUDE_DIRECTORIES>
  )

  # And then we link in the rest of the libraries just like in the shared case
  target_link_libraries(${target_name}
    PUBLIC
      Scine::UtilsOS
      outcome
    PRIVATE
      Boost::filesystem
      Boost::system
      RingDecomposerLibStatic
      nauty
      $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
  )

  add_eigen(${target_name} PUBLIC)
endfunction()

# Main library, of type as determined by BUILD_SHARED_LIBS
add_library(Molassembler
  $<TARGET_OBJECTS:molassembler_obj>
  ${CMAKE_CURRENT_BINARY_DIR}/Molassembler/Version.h
  ${CMAKE_CURRENT_BINARY_DIR}/Molassembler/Export.h
)
molassembler_includes_and_properties(Molassembler)
molassembler_library_links(Molassembler)
add_library(Scine::Molassembler ALIAS Molassembler)

if(BUILD_SHARED_LIBS)
  # If molassembler is shared, we still need a static variant for linking into
  # the tests and other binaries that hook deeper into molassembler than the
  # public API
  add_library(MolassemblerStatic STATIC
    $<TARGET_OBJECTS:molassembler_obj>
    ${CMAKE_CURRENT_BINARY_DIR}/Molassembler/Version.h
    ${CMAKE_CURRENT_BINARY_DIR}/Molassembler/Export.h
  )
  molassembler_includes_and_properties(MolassemblerStatic)
  molassembler_library_links(MolassemblerStatic)
  add_library(Scine::MolassemblerStatic ALIAS MolassemblerStatic)

  # Target name for others to hook into where needed
  set(MOLASSEMBLER_STATIC_TARGET MolassemblerStatic PARENT_SCOPE)
else()
  # If molassembler is static, then that's all we need
  set(MOLASSEMBLER_STATIC_TARGET Molassembler PARENT_SCOPE)
endif()

# Package version and config file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/ScineMolassemblerConfigVersion.cmake"
  VERSION ${Molassembler_VERSION}
  COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
  "config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/ScineMolassemblerConfig.cmake"
  INSTALL_DESTINATION "lib/cmake/ScineMolassembler"
)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/ScineMolassemblerConfigVersion.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/ScineMolassemblerConfig.cmake"
  DESTINATION "lib/cmake/ScineMolassembler"
)

# Target .cmake file
install(
  TARGETS Molassembler
  EXPORT molassemblerTargets
  DESTINATION lib
)

export(
  EXPORT molassemblerTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/ScineMolassemblerTargets.cmake"
  NAMESPACE Scine::
)

install(
  EXPORT molassemblerTargets
  FILE "ScineMolassemblerTargets.cmake"
  DESTINATION "lib/cmake/ScineMolassembler"
  NAMESPACE Scine::
)

# Headers
install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/Molassembler"
  DESTINATION include
  FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.hpp"
    PATTERN "*.hxx"
)

install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/Molassembler/Version.h
    ${CMAKE_CURRENT_BINARY_DIR}/Molassembler/Export.h
  DESTINATION include/Molassembler
)

# Source groups
source_group("Headers" FILES ${HEADERS_All})
source_group("Sources" FILES ${SOURCES_All})
