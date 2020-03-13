#
# This file is licensed under the 3-clause BSD license.
# Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
# See LICENSE.txt for details.
#
macro(import_utils_os)
  if(NOT TARGET SCINE::UtilsOS)
    # Try to find the package locally
    find_package(Scine OPTIONAL_COMPONENTS UtilsOS QUIET)
    if(NOT TARGET Scine::UtilsOS)
      cmessage(STATUS "Scine Utils not found, trying downloading instead...")
      # Download it instead
      include(DownloadProject)
      download_project(
        PROJ scine-utils-os
        GIT_REPOSITORY git@gitlab.chab.ethz.ch:scine/utils-open-source.git
        GIT_TAG        develop
        QUIET
      )

      # Note: Options defined in the project calling this function override default
      # option values specified in the imported project.
      set(_SCINE_BUILD_TESTS ${SCINE_BUILD_TESTS})
      set(SCINE_BUILD_TESTS OFF)
      add_subdirectory(${scine-utils-os_SOURCE_DIR} ${scine-utils-os_BINARY_DIR})
      set(SCINE_BUILD_TESTS ${_SCINE_BUILD_TESTS})
      unset(_SCINE_BUILD_TESTS)

      if(SCINE_BUILD_PYTHON_BINDINGS)
        set(SCINE_UTILS_PYTHON_BINARY_PATH ${scine-utils-os_BINARY_DIR}/src/Utils)
      endif()

      # Final check if all went well
      if(NOT TARGET Scine::UtilsOS)
        cmessage(FATAL_ERROR
          "Scine::UtilsOS was not found in your PATH and could not be "
          "established through a download. Try specifying Scine_DIR or "
          "altering CMAKE_PREFIX_PATH to point to a candidate Scine "
          "installation base directory."
        )
      endif()
    else()
      cmessage(STATUS "Scine Utilities found locally at ${Scine_DIR}")
    endif()
  endif()
endmacro()
