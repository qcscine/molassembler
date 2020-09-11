macro(import_rdl)
  # Do nothing if the target already exists
  if(TARGET RingDecomposerLib)
    return()
  endif()

  # Try to find as a package
  find_package(RingDecomposerLib QUIET)
  if(TARGET RingDecomposerLib)
    message(STATUS "RingDecomposerLib found locally at ${RingDecomposerLib_DIR}")
    return()
  endif()

  # Current state of the existing CMake code is unsuitable for subdirectory
  # inclusion and lacks installed config files, so prefer intrusive modeling
  include(DownloadProject)
  download_project(
    PROJ RingDecomposerLib
    GIT_REPOSITORY https://github.com/rareylab/RingDecomposerLib
    GIT_TAG v1.1.2
    QUIET
  )

  set(RDL_SOURCES_DIR "${RingDecomposerLib_SOURCE_DIR}/src/RingDecomposerLib")
  file(GLOB RDL_SOURCES "${RDL_SOURCES_DIR}/*.c")
  add_library(RingDecomposerLib STATIC ${RDL_SOURCES})
  target_include_directories(RingDecomposerLib PUBLIC
    $<BUILD_INTERFACE:${RDL_SOURCES_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
  )

  install(
    TARGETS RingDecomposerLib
    EXPORT RingDecomposerLibTargets
    ARCHIVE DESTINATION lib
  )
  install(FILES "${RDL_SOURCES_DIR}/RingDecomposerLib.h" DESTINATION include)

  include(CMakePackageConfigHelpers)
  write_basic_package_version_file(
    "${RingDecomposerLib_BINARY_DIR}/RingDecomposerLibConfigVersion.cmake"
    VERSION 1.1.3
    COMPATIBILITY AnyNewerVersion
  )

  file(WRITE "${RingDecomposerLib_BINARY_DIR}/config.cmake.in"
    "include(\$\{CMAKE_CURRENT_LIST_DIR\}/RingDecomposerLibTargets.cmake)\n @PACKAGE_INIT@"
  )

  configure_package_config_file(
    "${RingDecomposerLib_BINARY_DIR}/config.cmake.in"
    "${RingDecomposerLib_BINARY_DIR}/RingDecomposerLibConfig.cmake"
    INSTALL_DESTINATION "lib/cmake/RingDecomposerLib"
  )

  install(
    FILES
      "${RingDecomposerLib_BINARY_DIR}/RingDecomposerLibConfigVersion.cmake"
      "${RingDecomposerLib_BINARY_DIR}/RingDecomposerLibConfig.cmake"
    DESTINATION "lib/cmake/RingDecomposerLib"
  )

  export(
    EXPORT RingDecomposerLibTargets
    FILE "${RingDecomposerLib_BINARY_DIR}/RingDecomposerLibTargets.cmake"
  )

  install(
    EXPORT RingDecomposerLibTargets
    FILE "RingDecomposerLibTargets.cmake"
    DESTINATION "lib/cmake/RingDecomposerLib"
  )
endmacro()
