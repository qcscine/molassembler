#
# This file is licensed under the 3-clause BSD license.
# Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
# See LICENSE.txt for details.
#

macro(determine_lib_type target)
  get_target_property(_target_type ${target} TYPE)
  if(${_target_type} STREQUAL "SHARED_LIBRARY")
    set(LIB_TYPE "SHARED_LIBRARY")
  elseif(${_target_type} STREQUAL "STATIC_LIBRARY")
    set(LIB_TYPE "STATIC_LIBRARY")
  elseif(${_target_type} STREQUAL "UNKNOWN_LIBRARY")
    get_target_property(_target_location ${target} LOCATION)
    get_filename_component(_target_extension ${_target_location} EXT)
    if(${_target_extension} STREQUAL ${CMAKE_SHARED_LIBRARY_SUFFIX})
      set(LIB_TYPE "SHARED_LIBRARY")
    elseif(${_target_extension} STREQUAL ${CMAKE_STATIC_LIBRARY_SUFFIX})
      set(LIB_TYPE "STATIC_LIBRARY")
    else()
      set(LIB_TYPE "UNKNOWN_LIBRARY")
    endif()
    unset(_target_location)
    unset(_target_extension)
  endif()
  unset(_target_type)
endmacro()
