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

source_group("Headers" FILES ${HEADERS_All})

install(
  FILES LICENSE.txt
  DESTINATION share/doc/molassembler/licenses
  RENAME json.txt
)
