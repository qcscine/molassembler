# oldest tested version of cmake used (too lazy to figure out newest compatible
# cmake version)
cmake_minimum_required(VERSION 3.5.1)

# little variables
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# declaration
project(nlohmann_json_hpp)

set(json_VERSION 3.1.2)

set(HEADERS_All
  ${CMAKE_CURRENT_SOURCE_DIR}/json/json.hpp
)

add_library(json INTERFACE)
target_sources(json INTERFACE $<BUILD_INTERFACE:${HEADERS_All}>)
target_include_directories(
  json INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

install(
  TARGETS json
  EXPORT jsonTargets
  DESTINATION ${INCLUDE_INSTALL_DIR}
)

molassembler_install_component_cmake_files(
  COMPONENT json
  EXPORT_NAME jsonTargets
)

molassembler_install_component_headers(
  COMPONENT json
  DIRECTORY json
)

source_group("Headers" FILES ${HEADERS_All})
