#
# This file is licensed under the 3-clause BSD license.
# Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
# See LICENSE.txt for details.
#

macro(import_json)
  set(JSON_VERSION "3.9.1")
  set(JSON_LICENSE_FILE "${CMAKE_CURRENT_BINARY_DIR}/nlohmann/LICENSE")
  set(JSON_HPP_FILE "${CMAKE_CURRENT_BINARY_DIR}/nlohmann/json.hpp")

  include(DownloadFileHelper)
  try_resource_dir(
    SOURCE nlohmann/LICENSE nlohmann/json.hpp
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  )
  download_file(
    "https://raw.githubusercontent.com/nlohmann/json/v${JSON_VERSION}/single_include/nlohmann/json.hpp"
    ${JSON_HPP_FILE}
  )
  download_file(
    "https://raw.githubusercontent.com/nlohmann/json/v${JSON_VERSION}/LICENSE.MIT"
    ${JSON_LICENSE_FILE}
  )

  add_library(json INTERFACE)
  target_include_directories(json INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)

  install(
    FILES ${JSON_LICENSE_FILE}
    DESTINATION share/doc/molassembler/licenses
    RENAME nlohmann-json.txt
  )
endmacro()
