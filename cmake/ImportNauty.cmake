macro(import_nauty)
  # If it already exists as a target do nothing
  if(NOT TARGET nauty)
    # Try to find nauty as a package (for conan modeling since it doesn't have
    # its own CMake packaging files)
    find_package(nauty QUIET)
    if(TARGET nauty)
      message(STATUS "nauty found locally at ${nauty_DIR}")
    else()
      # Pull the archive if we haven't done so yet
      if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/nauty)
        include(DownloadFileHelper)
        download_file(
          "http://pallini.di.uniroma1.it/nauty27r1.tar.gz"
          ${CMAKE_CURRENT_BINARY_DIR}/nauty.tar.gz
        )
        # Unpack the archive and remove it
        execute_process(
          COMMAND ${CMAKE_COMMAND} -E tar zxf
          ${CMAKE_CURRENT_BINARY_DIR}/nauty.tar.gz
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        )
        file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/nauty.tar.gz)
        file(RENAME
          ${CMAKE_CURRENT_BINARY_DIR}/nauty27r1
          ${CMAKE_CURRENT_BINARY_DIR}/nauty
        )
      endif()

      set(NAUTY_HEADERS
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/gtools.h
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/gutils.h
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/naugroup.h
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/naugstrings.h
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/naurng.h
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/nausparse.h
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/nautaux.h
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/nautinv.h
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/naututil.h
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/nauty.h
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/planarity.h
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/quarticirred28.h
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/rng.h
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/schreier.h
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/traces.h
      )

      set(NAUTY_SOURCES
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/gtnauty.c
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/gtools.c
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/gutil1.c
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/gutil2.c
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/naugraph.c
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/naugroup.c
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/naurng.c
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/nausparse.c
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/nautil.c
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/nautinv.c
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/naututil.c
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/nauty.c
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/schreier.c
        ${CMAKE_CURRENT_BINARY_DIR}/nauty/traces.c
      )

      add_library(nauty ${NAUTY_HEADERS} ${NAUTY_SOURCES})
      target_include_directories(nauty PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
      )
      install(
        DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/nauty/
        DESTINATION include/nauty
        FILES_MATCHING PATTERN "*.h"
      )

      install(
        FILES ${CMAKE_CURRENT_BINARY_DIR}/nauty/COPYRIGHT
        DESTINATION share/doc/nauty/licenses
      )
      install(TARGETS nauty EXPORT nautyTargets ARCHIVE DESTINATION lib)

      include(CMakePackageConfigHelpers)
      write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/nauty/nautyConfigVersion.cmake"
        VERSION 2.7
        COMPATIBILITY AnyNewerVersion
      )

      file(
        WRITE "${CMAKE_CURRENT_BINARY_DIR}/nauty/config.cmake.in"
        "include(\$\{CMAKE_CURRENT_SOURCE_DIR\}/nautyTargets.cmake)\n@PACKAGE_INIT@"
      )

      configure_package_config_file(
        "${CMAKE_CURRENT_BINARY_DIR}/nauty/config.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/nauty/nautyConfig.cmake"
        INSTALL_DESTINATION "lib/cmake/nauty"
      )

      install(
        FILES
          "${CMAKE_CURRENT_BINARY_DIR}/nauty/nautyConfigVersion.cmake"
          "${CMAKE_CURRENT_BINARY_DIR}/nauty/nautyConfig.cmake"
        DESTINATION "lib/cmake/nauty"
      )

      export(
        EXPORT nautyTargets
        FILE "${CMAKE_CURRENT_BINARY_DIR}/nauty/nautyTargets.cmake"
      )

      install(
        EXPORT nautyTargets
        FILE "nautyTargets.cmake"
        DESTINATION "lib/cmake/nauty"
      )
    endif()
  endif()
endmacro()
