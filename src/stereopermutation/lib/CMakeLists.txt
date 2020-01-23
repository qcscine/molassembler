# dependencies
find_package(shapes REQUIRED)

set(HEADERS_All
  ${CMAKE_CURRENT_SOURCE_DIR}/stereopermutation/Composites.h
  ${CMAKE_CURRENT_SOURCE_DIR}/stereopermutation/Manipulation.h
  ${CMAKE_CURRENT_SOURCE_DIR}/stereopermutation/RotationEnumerator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/stereopermutation/Stereopermutation.h
)

set(SOURCES_All
  ${CMAKE_CURRENT_SOURCE_DIR}/stereopermutation/Composites.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/stereopermutation/Manipulation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/stereopermutation/RotationEnumerator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/stereopermutation/Stereopermutation.cpp
)

add_library(stereopermutation STATIC ${SOURCES_All} ${HEADERS_all})
target_include_directories(stereopermutation
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

target_include_directories(stereopermutation SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIRS})

target_link_libraries(stereopermutation PUBLIC shapes Eigen3::Eigen)

target_compile_options(stereopermutation PRIVATE ${MOLASSEMBLER_CXX_FLAGS})

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
