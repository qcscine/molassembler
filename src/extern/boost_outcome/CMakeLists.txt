# Included in source tree awaiting boost::outcome publication in Boost 1.69

# declaration
project(OutcomeProject)

set(boost_outcome_VERSION 0.1.0)

set(HEADERS_All
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/detail/outcome_exception_observers.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/detail/outcome_exception_observers_impl.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/detail/outcome_failure_observers.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/detail/result_error_observers.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/detail/result_final.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/detail/result_storage.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/detail/result_value_observers.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/detail/value_storage.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/policy/detail/common.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/policy/all_narrow.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/policy/outcome_error_code_throw_as_system_error.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/policy/outcome_exception_ptr_rethrow.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/policy/result_error_code_throw_as_system_error.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/policy/result_exception_ptr_rethrow.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/policy/terminate.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/policy/throw_bad_result_access.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/result.h
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/bad_access.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/config.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/convert.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/iostream_support.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/outcome.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/result.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/revision.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/success_failure.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/try.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/utils.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome/version.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost_outcome/outcome.hpp
)

add_library(boost_outcome INTERFACE)
target_sources(boost_outcome INTERFACE $<BUILD_INTERFACE:${HEADERS_All}>)
target_include_directories(
  boost_outcome INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
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
  DIRECTORY boost_outcome
)

source_group("Headers" FILES ${HEADERS_All})
