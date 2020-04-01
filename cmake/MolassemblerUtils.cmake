include(CMakePackageConfigHelpers)

set(INCLUDE_INSTALL_DIR "include")

function(molassembler_install_component_headers)
  set(bools "")
  set(valueArgs DIRECTORY COMPONENT)
  set(listArgs "")

  cmake_parse_arguments(ARGS "${bools}" "${valueArgs}" "${listArgs}" ${ARGN})

  install(
    DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${ARGS_DIRECTORY}"
    DESTINATION ${INCLUDE_INSTALL_DIR}
    COMPONENT ${ARGS_COMPONENT}
    FILES_MATCHING
      PATTERN "*.h"
      PATTERN "*.hpp"
      PATTERN "*.hxx"
  )
endfunction()

function(molassembler_install_component_cmake_files)
  set(bools "")
  set(valueArgs EXPORT_NAME COMPONENT)
  set(listArgs "")
  cmake_parse_arguments(ARGS "${bools}" "${valueArgs}" "${listArgs}" ${ARGN})

  # -- Config version file
  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${ARGS_COMPONENT}/${ARGS_COMPONENT}ConfigVersion.cmake"
    VERSION ${${ARGS_COMPONENT}_VERSION}
    COMPATIBILITY AnyNewerVersion
  )

  # -- Config file
  set(MOLASSEMBLER_CMAKE_ROOT "lib/cmake/molassembler")
  configure_package_config_file(
    "config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${ARGS_COMPONENT}/${ARGS_COMPONENT}Config.cmake"
    INSTALL_DESTINATION "${MOLASSEMBLER_CMAKE_ROOT}/${ARGS_COMPONENT}"
    PATH_VARS
      MOLASSEMBLER_CMAKE_ROOT
      INCLUDE_INSTALL_DIR
  )

  # Install config and configVersion
  install(
    FILES
      "${CMAKE_CURRENT_BINARY_DIR}/${ARGS_COMPONENT}/${ARGS_COMPONENT}Config.cmake"
      "${CMAKE_CURRENT_BINARY_DIR}/${ARGS_COMPONENT}/${ARGS_COMPONENT}ConfigVersion.cmake"
    DESTINATION "${MOLASSEMBLER_CMAKE_ROOT}/${ARGS_COMPONENT}"
    COMPONENT ${ARGS_COMPONENT}
  )

  # -- Targets
  # This makes the project targets importable from the build directory
  export(
    EXPORT ${ARGS_EXPORT_NAME}
    FILE "${CMAKE_CURRENT_BINARY_DIR}/${ARGS_COMPONENT}/${ARGS_COMPONENT}Targets.cmake"
    NAMESPACE molassembler::
  )

  # This makes the project targets importable from the install directory
  install(
    EXPORT ${ARGS_EXPORT_NAME}
    FILE "${ARGS_COMPONENT}Targets.cmake"
    DESTINATION "${MOLASSEMBLER_CMAKE_ROOT}/${ARGS_COMPONENT}"
    NAMESPACE molassembler::
    COMPONENT ${ARGS_COMPONENT}
  )

endfunction()
