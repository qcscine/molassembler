# Included in source tree awaiting boost::outcome publication in Boost 1.69

# oldest tested version of cmake used (too lazy to figure out newest compatible
# cmake version)
cmake_minimum_required(VERSION 3.5.1)

# little variables
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# declaration
project(OutcomeProject)

set(OUTCOME_VERSION 0.1.0)

set(HEADERS_All
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/detail/outcome_exception_observers.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/detail/outcome_exception_observers_impl.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/detail/outcome_failure_observers.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/detail/result_error_observers.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/detail/result_final.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/detail/result_storage.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/detail/result_value_observers.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/detail/value_storage.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/policy/detail/common.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/policy/all_narrow.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/policy/outcome_error_code_throw_as_system_error.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/policy/outcome_exception_ptr_rethrow.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/policy/result_error_code_throw_as_system_error.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/policy/result_exception_ptr_rethrow.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/policy/terminate.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/policy/throw_bad_result_access.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/result.h
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/bad_access.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/config.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/convert.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/iostream_support.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/outcome.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/result.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/revision.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/success_failure.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/try.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/utils.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome/version.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome/outcome.hpp
)

add_library(outcome INTERFACE)
target_sources(outcome INTERFACE $<BUILD_INTERFACE:${HEADERS_All}>)
target_include_directories(
  outcome INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

set(ConfigPackageLocation lib/cmake/boost-outcome)
set(IncludeInstallDir include/)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/boost-outcomeConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/boost-outcome/boost-outcomeConfig.cmake"
  INSTALL_DESTINATION ${ConfigPackageLocation}
  PATH_VARS IncludeInstallDir
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/boost-outcome/boost-outcomeConfigVersion.cmake"
  VERSION ${OUTCOME_VERSION}
  COMPATIBILITY AnyNewerVersion
)

install(
  TARGETS outcome
  EXPORT outcomeTargets
  DESTINATION include
)

install(
  EXPORT outcomeTargets
  DESTINATION "${ConfigPackageLocation}"
)

install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/boost-outcome"
  DESTINATION IncludeInstallDir
)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/boost-outcome/boost-outcomeConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/boost-outcome/boost-outcomeConfigVersion.cmake"
  DESTINATION ${ConfigPackageLocation}
)

source_group("Headers" FILES ${HEADERS_All})
