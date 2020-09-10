macro(import_outcome)
  include(DownloadFileHelper)

  set(OUTCOME_HEADER ${CMAKE_CURRENT_BINARY_DIR}/outcome/outcome.hpp)
  set(OUTCOME_VERSION "2.1.4")
  set(OUTCOME_LICENSE_FILE ${CMAKE_CURRENT_BINARY_DIR}/outcome/LICENSE)

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
    "include(\$\{CMAKE_CURRENT_SOURCE_DIR\}/OutcomeTargets.cmake)"
  )

  configure_package_config_file(
    "${CMAKE_CURRENT_BINARY_DIR}/outcome/config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/outcome/OutcomeConfig.cmake"
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
endmacro(import_outcome)
