cmake_minimum_required(VERSION 3.9)
project(nauty)

set(NAUTY_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/gtools.h
  ${CMAKE_CURRENT_SOURCE_DIR}/gutils.h
  ${CMAKE_CURRENT_SOURCE_DIR}/naugroup.h
  ${CMAKE_CURRENT_SOURCE_DIR}/naugstrings.h
  ${CMAKE_CURRENT_SOURCE_DIR}/naurng.h
  ${CMAKE_CURRENT_SOURCE_DIR}/nausparse.h
  ${CMAKE_CURRENT_SOURCE_DIR}/nautaux.h
  ${CMAKE_CURRENT_SOURCE_DIR}/nautinv.h
  ${CMAKE_CURRENT_SOURCE_DIR}/naututil.h
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty.h
  ${CMAKE_CURRENT_SOURCE_DIR}/planarity.h
  ${CMAKE_CURRENT_SOURCE_DIR}/quarticirred28.h
  ${CMAKE_CURRENT_SOURCE_DIR}/rng.h
  ${CMAKE_CURRENT_SOURCE_DIR}/schreier.h
  ${CMAKE_CURRENT_SOURCE_DIR}/traces.h
  )

set(NAUTY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/gtnauty.c
  ${CMAKE_CURRENT_SOURCE_DIR}/gtools.c
  ${CMAKE_CURRENT_SOURCE_DIR}/gutil1.c
  ${CMAKE_CURRENT_SOURCE_DIR}/gutil2.c
  ${CMAKE_CURRENT_SOURCE_DIR}/naugraph.c
  ${CMAKE_CURRENT_SOURCE_DIR}/naugroup.c
  ${CMAKE_CURRENT_SOURCE_DIR}/naurng.c
  ${CMAKE_CURRENT_SOURCE_DIR}/nausparse.c
  ${CMAKE_CURRENT_SOURCE_DIR}/nautil.c
  ${CMAKE_CURRENT_SOURCE_DIR}/nautinv.c
  ${CMAKE_CURRENT_SOURCE_DIR}/naututil.c
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty.c
  ${CMAKE_CURRENT_SOURCE_DIR}/schreier.c
  ${CMAKE_CURRENT_SOURCE_DIR}/traces.c
  )

add_library(nauty STATIC ${NAUTY_HEADERS} ${NAUTY_SOURCES})
target_include_directories(nauty PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
  )
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  DESTINATION include/nauty
  FILES_MATCHING PATTERN "*.h"
  )

install(
  FILES ${CMAKE_CURRENT_SOURCE_DIR}/COPYRIGHT
  DESTINATION share/doc/nauty/licenses
  )
install(TARGETS nauty EXPORT nautyTargets ARCHIVE DESTINATION lib)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/nauty/nautyConfigVersion.cmake"
  VERSION 2.7
  COMPATIBILITY AnyNewerVersion
  )

file(
  WRITE "${CMAKE_CURRENT_BINARY_DIR}/nauty/config.cmake.in"
  "include(\$\{CMAKE_CURRENT_LIST_DIR\}/nautyTargets.cmake)\n@PACKAGE_INIT@"
  )

configure_package_config_file(
  "${CMAKE_CURRENT_BINARY_DIR}/nauty/config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/nauty/nautyConfig.cmake"
  INSTALL_DESTINATION "lib/cmake/nauty"
  )

install(
  FILES
  "${CMAKE_CURRENT_BINARY_DIR}/nauty/nautyConfigVersion.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/nauty/nautyConfig.cmake"
  DESTINATION "lib/cmake/nauty"
  )

export(
  EXPORT nautyTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/nauty/nautyTargets.cmake"
  )

install(
  EXPORT nautyTargets
  FILE "nautyTargets.cmake"
  DESTINATION "lib/cmake/nauty"
  )
