project(nauty VERSION 2.6)

set(NAUTY_HEADER
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/gtools.h
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/gutils.h
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/naugroup.h
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/naugstrings.h
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/naurng.h
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/nausparse.h
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/nautaux.h
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/nautinv.h
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/naututil.h
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/nauty.h
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/planarity.h
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/quarticirred28.h
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/rng.h
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/schreier.h
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/traces.h
)

set(NAUTY_CPP
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/gtnauty.c
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/gtools.c
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/gutil1.c
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/gutil2.c
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/naugraph.c
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/naugroup.c
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/naurng.c
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/nausparse.c
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/nautil.c
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/nautinv.c
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/naututil.c
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/nauty.c
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/schreier.c
  ${CMAKE_CURRENT_SOURCE_DIR}/nauty/traces.c
)

add_library(nauty STATIC ${NAUTY_HEADER} ${NAUTY_CPP})
target_include_directories(nauty PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)
set_target_properties(nauty PROPERTIES POSITION_INDEPENDENT_CODE ON)
install(TARGETS nauty EXPORT NautyTargets ARCHIVE DESTINATION lib)

molassembler_install_component_cmake_files(
  COMPONENT nauty
  EXPORT_NAME NautyTargets
)
