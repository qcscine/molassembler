#
# This file is licensed under the 3-clause BSD license.
# Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
# See LICENSE.txt for details.
#

macro(import_rdl)
  # Do nothing if the target already exists
  if(NOT TARGET RingDecomposerLib)
    # Try to find as a package
    find_package(RingDecomposerLib QUIET)
    if(TARGET RingDecomposerLib)
      message(STATUS "RingDecomposerLib found locally at ${RingDecomposerLib_DIR}")
    else()
      include(DownloadFileHelper)
      if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/RingDecomposerLib)
        try_resource_dir(
          SOURCE v1.1.2.tar.gz
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
        )
        download_file(
          "https://github.com/rareylab/RingDecomposerLib/archive/v1.1.2.tar.gz"
          ${CMAKE_CURRENT_BINARY_DIR}/v1.1.2.tar.gz
        )
        execute_process(
          COMMAND ${CMAKE_COMMAND} -E tar zxf
          ${CMAKE_CURRENT_BINARY_DIR}/v1.1.2.tar.gz
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        )
        file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/v1.1.2.tar.gz)
        file(RENAME
          ${CMAKE_CURRENT_BINARY_DIR}/RingDecomposerLib-1.1.2
          ${CMAKE_CURRENT_BINARY_DIR}/RingDecomposerLib
        )
      endif()

      # Current state of the existing CMake code is unsuitable for subdirectory
      # inclusion and lacks installed config files, so prefer intrusive modeling
      set(RDL_ROOT ${CMAKE_CURRENT_BINARY_DIR}/RingDecomposerLib)
      set(RDL_SOURCES_DIR "${RDL_ROOT}/src/RingDecomposerLib")
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
        "${RDL_ROOT}/RingDecomposerLibConfigVersion.cmake"
        VERSION 1.1.2
        COMPATIBILITY AnyNewerVersion
      )

    file(WRITE "${RDL_ROOT}/config.cmake.in"
        "include(\$\{CMAKE_CURRENT_LIST_DIR\}/RingDecomposerLibTargets.cmake)\n @PACKAGE_INIT@"
      )

      configure_package_config_file(
        "${RDL_ROOT}/config.cmake.in"
        "${RDL_ROOT}/RingDecomposerLibConfig.cmake"
        INSTALL_DESTINATION "lib/cmake/RingDecomposerLib"
      )

      install(
        FILES
          "${RDL_ROOT}/RingDecomposerLibConfigVersion.cmake"
          "${RDL_ROOT}/RingDecomposerLibConfig.cmake"
        DESTINATION "lib/cmake/RingDecomposerLib"
      )

      export(
        EXPORT RingDecomposerLibTargets
        FILE "${RDL_ROOT}/RingDecomposerLibTargets.cmake"
      )

      install(
        EXPORT RingDecomposerLibTargets
        FILE "RingDecomposerLibTargets.cmake"
        DESTINATION "lib/cmake/RingDecomposerLib"
      )
    endif()
  endif()
endmacro()
