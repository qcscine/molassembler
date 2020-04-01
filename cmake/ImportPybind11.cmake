#
# This file is licensed under the 3-clause BSD license.
# Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
# See LICENSE.txt for details.
#
macro(import_pybind11)
  # If the target already exists, do nothing
  if(NOT TARGET pybind11::pybind11)
    find_package(pybind11 2.4.2 EXACT QUIET)
    if(TARGET pybind11::pybind11)
      cmessage(STATUS "Found pybind11 at ${pybind11_DIR}")
    else()
      cmessage(STATUS "Pybind11 was not found, trying to download it instead...")
      # Download it instead
      include(DownloadProject)
      download_project(
        PROJ pybind11
        GIT_REPOSITORY      https://github.com/pybind/pybind11.git
        GIT_TAG             v2.4.2
        QUIET
      )

      add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})

      # Final check if all went well
      if(NOT EXISTS "${pybind11_SOURCE_DIR}/CMakeLists.txt")
        cmessage(FATAL_ERROR "Pybind11 could not be established through a download.")
      endif()
    endif()
  endif()
endmacro()
