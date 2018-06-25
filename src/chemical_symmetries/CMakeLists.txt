# declaration
project(chemicalsymmetries VERSION 0.1.0)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # Flags
  # Clang constexpr options
  set(CONSTEXPR_OPTIONS "-fconstexpr-steps=1999999999 -fconstexpr-backtrace-limit=0")
else()
  string(
    CONCAT buildWarning
    "Due to compiler restrictions, this build of chemical_symmetries will not"
    " include the precomputation of constexpr properties. This will incur "
    "a runtime speed penalty. Consider building it with Clang instead."
  )
  message(WARNING ${buildWarning})
endif()

enable_testing()

add_subdirectory(lib)
add_subdirectory(test)
add_subdirectory(analysis)
