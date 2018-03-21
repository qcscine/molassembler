cmake_minimum_required(VERSION 3.5)

project(molassembler)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  string(
    CONCAT FAILURE_MESSAGE
    "molassembler cannot be built with GCC. Building"
    " the chemical_symmetries library fails, perhaps because of GCC "
    "bug 71504: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=71504\n"
    "Clang builds work fine!"
  )
  message(FATAL_ERROR ${FAILURE_MESSAGE})
else()
  enable_testing()
  add_subdirectory(src)
endif()
