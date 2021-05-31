#
# This file is licensed under the 3-clause BSD license.
# Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
# See LICENSE.txt for details.
#

# Use C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Generate a compile-commands database for use with linters / clang-tidy
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Generate position independent code in all targets
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Handle opportunity for interprocedural optimization
if(MOLASSEMBLER_IPO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT canLTO OUTPUT LTOOut)
  if(canLTO)
    # Set default value of IPO for all targets produced
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
    cmessage(STATUS "Interprocedural optimization has been enabled")
  else()
    cmessage(WARNING "Interprocedural optimization is not supported: ${LTOOut}")
  endif()
endif()

# Handle opportunity for architecture-specific instruction sets
if(NOT "${SCINE_MARCH}" STREQUAL "" AND NOT MSVC)
  list(APPEND MOLASSEMBLER_CXX_FLAGS -march=${SCINE_MARCH})
endif()

# Compilation flags
if(MSVC)
  # Compilation on different cores (faster)
  list(APPEND MOLASSEMBLER_CXX_FLAGS /MP)
  # We require way more variadic templates in the constexpr evaluation
  # No idea if this is still relevant, but we're going to assume so
  list(APPEND MOLASSEMBLER_CXX_FLAGS /D _VARIADIC_MAX=200)
  # Allow more constexpr steps
  list(APPEND MOLASSEMBLER_CXX_FLAGS /constexpr:steps 100000000) # 100M
  # Enable use of math constants like M_PI, M_LN10 by compiling sources with
  add_definitions(/D_USE_MATH_DEFINES)
else()
  # This ought to work for GCC and Clang equally
  list(APPEND MOLASSEMBLER_CXX_FLAGS -Wall -Wpedantic -Wextra -Wshadow)

  # Some GCC-specific compiler options
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    cmessage(STATUS "Enabling GCC specific warning flags")
    # Handle no gnu unique flags
    if(MOLASSEMBLER_NO_GNU_UNIQUE)
      list(APPEND MOLASSEMBLER_CXX_FLAGS --no-gnu-unique)
    endif()

    list(APPEND MOLASSEMBLER_CXX_FLAGS
      -Wduplicated-cond
      -Wlogical-op
      -Wold-style-cast
      -Wuseless-cast
      -Wdouble-promotion
      -Wno-maybe-uninitialized
    )
  endif()

  if(
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"
    OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang"
  )
    list(APPEND MOLASSEMBLER_CXX_FLAGS
      -fconstexpr-steps=100000000
      -fconstexpr-backtrace-limit=0
    )
  endif()
endif()

# Force the linker to report undefined symbols in shared libraries
if(UNIX AND NOT APPLE)
  string(APPEND CMAKE_SHARED_LINKER_FLAGS " -Wl,--no-undefined")
endif()
