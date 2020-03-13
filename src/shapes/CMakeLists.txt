if(
  NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"
  OR NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang"
)
  cmessage(WARNING
    "Due to compiler restrictions, shapes will not "
    "precompute constexpr properties."
  )
endif()

add_subdirectory(lib)
add_subdirectory(test)
add_subdirectory(analysis)
