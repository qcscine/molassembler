set(stereopermutation_VERSION 0.1.0)

# dependencies
find_package(Eigen3 3.1.1 REQUIRED)
find_package(chemical_symmetries REQUIRED)

set(HEADERS_All
  ${CMAKE_CURRENT_SOURCE_DIR}/stereopermutation/Stereopermutation.h
  ${CMAKE_CURRENT_SOURCE_DIR}/stereopermutation/GenerateUniques.h
  ${CMAKE_CURRENT_SOURCE_DIR}/stereopermutation/LogicalOperatorTests.h
  ${CMAKE_CURRENT_SOURCE_DIR}/stereopermutation/Util.h
)

set(SOURCES_All
  ${CMAKE_CURRENT_SOURCE_DIR}/stereopermutation/Stereopermutation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/stereopermutation/GenerateUniques.cpp
)

add_library(stereopermutation STATIC ${SOURCES_All} ${HEADERS_all})
target_include_directories(stereopermutation PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_include_directories(stereopermutation SYSTEM PUBLIC $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>)
target_include_directories(stereopermutation PUBLIC ${EIGEN3_INCLUDE_DIRS})

if(TARGET chemical_symmetries)
  message(STATUS "chemical_symmetries is a target")
else()
  message(STATUS "chemical_symmetries isn't a target")
endif()

target_link_libraries(stereopermutation PUBLIC chemical_symmetries)

if(UNIX)
  add_definitions(-fPIC)
endif()

target_compile_options(stereopermutation PUBLIC -Wall -Wpedantic -Wextra)

if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  message(STATUS "Selected release build")
  # no debug information, disable assert
  target_compile_definitions(stereopermutation PRIVATE NDEBUG)
  target_compile_options(stereopermutation PRIVATE -O2)
else() # assume Debug
  message(STATUS "Assume debug build")
  target_compile_options(stereopermutation PRIVATE -g -O2)
endif()

install(
  TARGETS stereopermutation
  EXPORT stereopermutationTargets
  DESTINATION lib
)

molassembler_install_component_cmake_files(
  COMPONENT stereopermutation
  EXPORT_NAME stereopermutationTargets
)

molassembler_install_component_headers(
  COMPONENT stereopermutation
  DIRECTORY stereopermutation
)

# Source groups for Visual Studio
source_group("Headers" FILES ${HEADERS_All})
source_group("Sources" FILES ${SOURCES_All})
