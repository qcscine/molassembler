# overwrite find_package with test against subproject list
set(as_subproject shapes cyclic_polygons stereopermutation temple)
macro(find_package)
  if(NOT "${ARGV0}" IN_LIST as_subproject)
    _find_package(${ARGV})
  endif()
endmacro()

# CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(ColorMessages)

if(MOLASSEMBLER_SANITIZE)
  cmessage(STATUS "MOLASSEMBLER_SANITIZE: Adding address and UB sanitizers")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address,undefined")
  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -fsanitize=address,undefined")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address,undefined")

  cmessage(WARNING
    "Address sanitiziation can cause ABI incompatibility issues "
    "if dependencies have not been compiled the same way."
  )
endif()

if(COVERAGE AND CMAKE_BUILD_TYPE STREQUAL "Debug")
  cmessage(STATUS "COVERAGE: Adding coverage flags")
  set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -coverage ${CMAKE_CXX_FLAGS_DEBUG}")
endif()

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
