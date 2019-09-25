project(molassembler)

# project declaration
enable_testing()

include(ImportUtilsOS)
import_utils_os()
add_subdirectory(lib)
add_subdirectory(test)
add_subdirectory(analysis)
add_subdirectory(validation)

if(SCINE_BUILD_PYTHON_BINDINGS)
  add_subdirectory(python)
endif()

# Copy data files to build directory for testing
file(
  COPY ${CMAKE_CURRENT_SOURCE_DIR}/data
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)
