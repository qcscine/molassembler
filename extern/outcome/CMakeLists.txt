# Included in source tree awaiting boost::outcome release in Boost 1.69

# declaration
project(Outcome VERSION 2.0.0)

set(OUTCOME_SINGLE_HEADER_FILE ${CMAKE_CURRENT_BINARY_DIR}/outcome/outcome.hpp)
if(NOT EXISTS ${OUTCOME_SINGLE_HEADER_FILE})
  cmessage(STATUS "Pulling boost outcome single header...")
  file(DOWNLOAD
    https://github.com/ned14/outcome/raw/master/single-header/outcome.hpp
    ${CMAKE_CURRENT_BINARY_DIR}/outcome/outcome.hpp
    STATUS OUTCOME_PULL_STATUS
  )
  list(GET OUTCOME_PULL_STATUS 0 OUTCOME_PULL_RESULT)
  if(NOT ${OUTCOME_PULL_RESULT} EQUAL 0)
    if(EXISTS ${OUTCOME_SINGLE_HEADER_FILE})
      file(REMOVE ${OUTCOME_SINGLE_HEADER_FILE})
    endif()
    list(GET OUTCOME_PULL_STATUS 1 OUTCOME_PULL_ERROR)
    cmessage(FATAL_ERROR
      "Could not pull outcome header: ${OUTCOME_PULL_ERROR}"
      "Please download the header at "
      "https://github.com/ned14/outcome/raw/master/single-header/outcome.hpp "
      "and place it in your build tree at "
      "build/src/extern/outcome/outcome/outcome.hpp"
    )
  endif()
endif()

set(OUTCOME_HEADER ${CMAKE_CURRENT_BINARY_DIR}/outcome/outcome.hpp)

add_library(outcome INTERFACE)
target_sources(outcome INTERFACE $<BUILD_INTERFACE:${OUTCOME_HEADER}>)
target_include_directories(outcome INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

install(TARGETS outcome EXPORT outcomeTargets)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/OutcomeConfigVersion.cmake"
  VERSION ${Outcome_VERSION}
  COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
  "config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/OutcomeConfig.cmake"
  INSTALL_DESTINATION "lib/cmake/Outcome"
)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/OutcomeConfigVersion.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/OutcomeConfig.cmake"
  DESTINATION "lib/cmake/Outcome"
)

export(
  EXPORT outcomeTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/OutcomeTargets.cmake"
)

install(
  EXPORT outcomeTargets
  FILE "OutcomeTargets.cmake"
  DESTINATION "lib/cmake/Outcome"
)

install(FILES ${OUTCOME_HEADER} DESTINATION include/outcome)

source_group("Headers" FILES ${OUTCOME_HEADER})
