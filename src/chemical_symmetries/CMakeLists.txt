# declaration
project(chemicalsymmetries VERSION 0.1.0)

if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  string(
    CONCAT buildWarning
    "Due to compiler restrictions, this build of chemical_symmetries will not"
    " include the precomputation of constexpr properties. This will incur "
    "a runtime speed penalty. Clang is required for the precomputation."
  )
  message(WARNING ${buildWarning})
endif()

enable_testing()

add_subdirectory(lib)
add_subdirectory(test)
add_subdirectory(analysis)
