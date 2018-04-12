# oldest tested version of cmake used (too lazy to figure out newest compatible
# cmake version)
cmake_minimum_required(VERSION 3.5.1)

# declaration
project(chemical_symmetries)

set(WARNING_FLAGS "-Wall -Wpedantic -Wextra")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # Flags
  # Clang constexpr options
  set(CONSTEXPR_OPTIONS "-fconstexpr-steps=1999999999 -fconstexpr-backtrace-limit=0")

  enable_testing()

  add_subdirectory(lib)
  add_subdirectory(test)
  add_subdirectory(analysis)
else()
  string(
    CONCAT buildWarning
    "Due to compiler restrictions, this build of chemical_symmetries will not"
    " include the precomputation of constexpr properties. This will incur "
    "a runtime speed penalty."
  )
  message(WARNING ${buildWarning})
endif()
