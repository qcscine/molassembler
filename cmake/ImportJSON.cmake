macro(import_json)
  include(DownloadFileHelper)

  set(JSON_VERSION "3.9.1")
  set(JSON_LICENSE_FILE "${CMAKE_CURRENT_BINARY_DIR}/nlohmann/LICENSE")

  download_file(
    "https://raw.githubusercontent.com/nlohmann/json/v${JSON_VERSION}/single_include/nlohmann/json.hpp"
    "${CMAKE_CURRENT_BINARY_DIR}/nlohmann/json.hpp"
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
