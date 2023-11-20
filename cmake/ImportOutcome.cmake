#
# This file is licensed under the 3-clause BSD license.
# Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
# See LICENSE.txt for details.
#

macro(import_outcome)
  if(NOT DEFINED Boost_VERSION)
    find_package(Boost QUIET)
  endif()
  
  set(BOOST_HAS_OUTCOME 0)
  if(${Boost_VERSION} MATCHES "\\.")
    if(${Boost_VERSION} VERSION_GREATER_EQUAL 1.70.0)
       set(BOOST_HAS_OUTCOME 1)
    endif()
  else()
    if(${Boost_VERSION} GREATER_EQUAL 107000)
       set(BOOST_HAS_OUTCOME 1)
    endif()
  endif()

  if(BOOST_HAS_OUTCOME)
    set(OUTCOME_NAMESPACE "BOOST_OUTCOME_V2_NAMESPACE")
    set(OUTCOME_INCLUDE_HEADER "boost/outcome.hpp")
  else()
    include(DownloadFileHelper)

    set(OUTCOME_HEADER ${CMAKE_CURRENT_BINARY_DIR}/outcome/outcome.hpp)
    set(OUTCOME_VERSION "2.2.4")
    set(OUTCOME_LICENSE_FILE ${CMAKE_CURRENT_BINARY_DIR}/outcome/LICENSE)

    try_resource_dir(
      SOURCE outcome/outcome.hpp outcome/LICENSE
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
    )

    download_file(
      "https://raw.githubusercontent.com/ned14/outcome/v${OUTCOME_VERSION}/single-header/outcome.hpp"
      ${OUTCOME_HEADER}
    )

    download_file(
      "https://raw.githubusercontent.com/ned14/outcome/v${OUTCOME_VERSION}/Licence.txt"
      ${OUTCOME_LICENSE_FILE}
    )

    add_library(outcome INTERFACE)
    target_sources(outcome INTERFACE $<BUILD_INTERFACE:${OUTCOME_HEADER}>)
    target_include_directories(outcome INTERFACE
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
      $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
    )

    install(TARGETS outcome EXPORT outcomeTargets)

    include(CMakePackageConfigHelpers)
    write_basic_package_version_file(
      "${CMAKE_CURRENT_BINARY_DIR}/outcome/OutcomeConfigVersion.cmake"
      VERSION ${OUTCOME_VERSION}
      COMPATIBILITY AnyNewerVersion
    )

    file(
      WRITE "${CMAKE_CURRENT_BINARY_DIR}/outcome/config.cmake.in"
      "include(\$\{CMAKE_CURRENT_LIST_DIR\}/OutcomeTargets.cmake)\n@PACKAGE_INIT@"
    )

    configure_package_config_file(
      "${CMAKE_CURRENT_BINARY_DIR}/outcome/config.cmake.in"
      "${CMAKE_CURRENT_BINARY_DIR}/outcome/OutcomeConfig.cmake"
      INSTALL_DESTINATION "lib/cmake/Outcome"
    )

    install(
      FILES
        "${CMAKE_CURRENT_BINARY_DIR}/outcome/OutcomeConfigVersion.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/outcome/OutcomeConfig.cmake"
      DESTINATION "lib/cmake/Outcome"
    )

    export(
      EXPORT outcomeTargets
      FILE "${CMAKE_CURRENT_BINARY_DIR}/outcome/OutcomeTargets.cmake"
    )

    install(
      EXPORT outcomeTargets
      FILE "OutcomeTargets.cmake"
      DESTINATION "lib/cmake/Outcome"
    )

    install(
      FILES ${OUTCOME_LICENSE_FILE}
      DESTINATION share/doc/molassembler/licenses
      RENAME outcome.txt
    )

    install(FILES ${OUTCOME_HEADER} DESTINATION include/outcome)

    set(OUTCOME_NAMESPACE "OUTCOME_V2_NAMESPACE")
    set(OUTCOME_INCLUDE_HEADER "outcome/outcome.hpp")
  endif()
endmacro(import_outcome)
