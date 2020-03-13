# Included in source tree awaiting boost::outcome publication in Boost 1.69

# declaration
project(Outcome)

set(boost_outcome_VERSION 2.0.0)

set(OUTCOME_SINGLE_HEADER_FILE ${CMAKE_CURRENT_BINARY_DIR}/boost_outcome/outcome.hpp)
if(NOT EXISTS ${OUTCOME_SINGLE_HEADER_FILE})
  cmessage(STATUS "Pulling boost outcome single header...")
  file(DOWNLOAD
    https://github.com/ned14/outcome/raw/master/single-header/outcome.hpp
    ${CMAKE_CURRENT_BINARY_DIR}/boost_outcome/outcome.hpp
  )
endif()

set(OUTCOME_HEADER ${CMAKE_CURRENT_BINARY_DIR}/boost_outcome/outcome.hpp)

add_library(boost_outcome INTERFACE)
target_sources(boost_outcome INTERFACE $<BUILD_INTERFACE:${OUTCOME_HEADER}>)
target_include_directories(boost_outcome INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

install(
  TARGETS boost_outcome
  EXPORT boost_outcomeTargets
  DESTINATION ${INCLUDE_INSTALL_DIR}
)

molassembler_install_component_cmake_files(
  COMPONENT boost_outcome
  EXPORT_NAME boost_outcomeTargets
)

molassembler_install_component_headers(
  COMPONENT boost_outcome
  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/boost_outcome
)

source_group("Headers" FILES ${OUTCOME_HEADER})
