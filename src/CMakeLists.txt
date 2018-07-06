set(molassembler_VERSION_MAJOR 0)
set(molassembler_VERSION_MINOR 1)
set(molassembler_VERSION_PATCH 0)

set(as_subproject chemical_symmetries cyclic_polygons stereopermutation temple)

# Project-specific options
option(MOLASSEMBLER_IPO "Try to enable interprocedural optimization" ON)

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
include(Utils)

# Permit limited builds
if(DEFINED MOLASSEMBLER_INCLUDE_PROJECT_DIRS)
  foreach(subdirectory ${MOLASSEMBLER_INCLUDE_PROJECT_DIRS})
    add_subdirectory(${subdirectory})
  endforeach()
else()
  add_subdirectory(extern)
  add_subdirectory(temple)
  add_subdirectory(cyclic_polygons)
  add_subdirectory(gor1)
  add_subdirectory(chemical_symmetries)
  add_subdirectory(stereopermutation)
  add_subdirectory(molassembler)
endif()

# Add molassembler global configuration
install(
  FILES "${CMAKE_CURRENT_SOURCE_DIR}/cmake/molassemblerConfig.cmake"
  DESTINATION lib/cmake/molassembler
)
