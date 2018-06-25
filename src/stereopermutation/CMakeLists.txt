project(stereopermutations)

enable_testing()

message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "Compiler version: ${CMAKE_CXX_COMPILER_VERSION}")

add_subdirectory(lib)
add_subdirectory(test)
add_subdirectory(analysis)
