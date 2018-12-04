function(scine_install_component_cmake_files)
  set(options "")
  set(oneValueArgs COMPONENT EXPORT_NAME CONFIG_FILE)
  set(multiValueArgs "")
  cmake_parse_arguments(SCINE_INSTALL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  include(CMakePackageConfigHelpers)

  # -- Config version file
  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${SCINE_INSTALL_COMPONENT}/${SCINE_INSTALL_COMPONENT}ConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
  )

  set(SCINE_COMPONENT ${PROJECT_NAME})
  set(SCINE_COMPONENT_VERSION ${PROJECT_VERSION})

  # -- Config file
  configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${SCINE_INSTALL_COMPONENT}/${SCINE_INSTALL_COMPONENT}Config.cmake"
    INSTALL_DESTINATION "${SCINE_CMAKE_PACKAGE_ROOT}/${SCINE_INSTALL_COMPONENT}"
    PATH_VARS SCINE_CMAKE_PACKAGE_ROOT
  )

  # Install config and configVersion
  install(
    FILES
      "${CMAKE_CURRENT_BINARY_DIR}/${SCINE_INSTALL_COMPONENT}/${SCINE_INSTALL_COMPONENT}Config.cmake"
      "${CMAKE_CURRENT_BINARY_DIR}/${SCINE_INSTALL_COMPONENT}/${SCINE_INSTALL_COMPONENT}ConfigVersion.cmake"
    DESTINATION "${SCINE_CMAKE_PACKAGE_ROOT}/${SCINE_INSTALL_COMPONENT}"
    COMPONENT ${SCINE_INSTALL_COMPONENT}
  )

  # -- Targets file
  # This makes the project importable from the build directory
  export(
    EXPORT ${SCINE_INSTALL_EXPORT_NAME}
    FILE "${CMAKE_CURRENT_BINARY_DIR}/${SCINE_INSTALL_COMPONENT}/${SCINE_INSTALL_COMPONENT}Targets.cmake"
    NAMESPACE Scine::
  )

  # This makes the project importable from the install directory
  # Put config file in per-project dir (name MUST match), can also
  # just go into <prefix>/cmake.
  install(
    EXPORT ${SCINE_INSTALL_EXPORT_NAME}
    FILE "${SCINE_INSTALL_COMPONENT}Targets.cmake"
    DESTINATION "${SCINE_CMAKE_PACKAGE_ROOT}/${SCINE_INSTALL_COMPONENT}"
    NAMESPACE Scine::
    COMPONENT ${SCINE_INSTALL_COMPONENT}
  )
endfunction()

function(workaround_link_object_library_target targetName)
  set(options "")
  set(oneValueArgs "")
  set(multiValueArgs PUBLIC PRIVATE)
  cmake_parse_arguments(WORKAROUND_LINK "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.12.0")
    target_link_libraries(targetName
      PUBLIC ${WORKAROUND_LINK_PUBLIC}
      PRIVATE ${WORKAROUND_LINK_PRIVATE}
    )
  else()
    foreach(publicTarget ${WORKAROUND_LINK_PUBLIC})
      target_include_directories(targetName
        PUBLIC $<TARGET_PROPERTY:${publicTarget},INTERFACE_INCLUDE_DIRECTORIES>
      )
    endforeach()
    foreach(privateTarget ${WORKAROUND_LINK_PRIVATE})
      target_include_directories(targetName
        PRIVATE $<TARGET_PROPERTY:${privateTarget},INTERFACE_INCLUDE_DIRECTORIES>
      )
    endforeach()
  endif()
endfunction()
