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
    message(STATUS "Interprocedural optimization has been enabled")
  else()
    message(WARNING "Interprocedural optimization is not supported: ${LTOOut}")
  endif()
endif()

# Compilation flags
if(MSVC)
  # Compilation on different cores (faster)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
  # We require way more variadic templates in the constexpr evaluation
  # No idea if this is still relevant, but we're going to assume so
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D _VARIADIC_MAX=200")
  # Allow more constexpr steps
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /constexpr:steps:10000000") # 10M
else()
  # This ought to work for GCC and Clang equally
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic -Wextra -Wshadow")

  # Some GCC-specific compiler options
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message(STATUS "Enabling GCC specific warning flags")
    # Handle no gnu unique flags
    if(MOLASSEMBLER_NO_GNU_UNIQUE)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --no-gnu-unique")
    endif()

    set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} \
      -Wsuggest-attribute=const \
      -Wsuggest-attribute=pure \
      -Wsuggest-final-types \
      -Wsuggest-final-methods"
    )
  endif()
endif()

# Force the linker to report undefined symbols in shared libraries
if(UNIX AND NOT APPLE)
  string(APPEND CMAKE_SHARED_LINKER_FLAGS " -Wl,--no-undefined")
endif()
