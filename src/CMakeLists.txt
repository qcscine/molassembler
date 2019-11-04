# overwrite find_package with test against subproject list
set(as_subproject shapes cyclic_polygons stereopermutation temple)
macro(find_package)
  if(NOT "${ARGV0}" IN_LIST as_subproject)
    _find_package(${ARGV})
  endif()
endmacro()

# CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(ImportUtilsOS)
import_utils_os()

# Externalize parts for sake of brevity
include(CompilerSpecificSetup)
include(Dependencies)
include(MolassemblerUtils)

# Subdirectories
add_subdirectory(extern)
add_subdirectory(temple)
add_subdirectory(cyclic_polygons)
add_subdirectory(gor1)
add_subdirectory(shapes)
add_subdirectory(stereopermutation)
add_subdirectory(molassembler)

# Add molassembler global configuration
install(
  FILES "${CMAKE_CURRENT_SOURCE_DIR}/cmake/molassemblerConfig.cmake"
  DESTINATION lib/cmake/molassembler
)
