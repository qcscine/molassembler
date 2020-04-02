# Check for Doxygen, if not found, return
find_package(Doxygen)
if(NOT DOXYGEN_FOUND)
  message(WARNING "Doxygen was not found. The C++ documentation cannot be built.")
  return()
endif()

# Set up the doxygen build
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
  ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  @ONLY
)

add_custom_target(doc
  ALL
  COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generate documentation with Doxygen" VERBATIM
)

install(
  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/
  DESTINATION share/doc/molassembler/cpp
  MESSAGE_NEVER
  OPTIONAL
)
