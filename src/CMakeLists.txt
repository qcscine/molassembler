set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(molassembler_VERSION_MAJOR 0)
set(molassembler_VERSION_MINOR 1)
set(molassembler_VERSION_PATCH 0)

set(as_subproject chemical_symmetries cyclic_polygons stereopermutation temple)

# overwrite find_package with test against subproject list
macro(find_package)
  if(NOT "${ARGV0}" IN_LIST as_subproject)
    _find_package(${ARGV})
  endif()
endmacro()

# CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Externalize parts for sake of brevity
include(CompilerSpecificSetup)
include(Dependencies)

add_subdirectory(temple)
add_subdirectory(cyclic_polygons)
add_subdirectory(gor1)
add_subdirectory(boost-outcome)
add_subdirectory(chemical_symmetries)
add_subdirectory(stereopermutation)
add_subdirectory(molassembler)
