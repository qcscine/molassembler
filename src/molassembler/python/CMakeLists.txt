include(ImportPybind11)
import_pybind11()
set(PYBIND11_PYTHON_VERSION ${PYTHONVERSION})

file(GLOB_RECURSE MOLASSEMBLER_PYTHON_CPPS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# Python module
pybind11_add_module(molassembler_python ${MOLASSEMBLER_PYTHON_CPPS})
set_target_properties(molassembler_python PROPERTIES
  SUFFIX ".so"
  LIBRARY_OUTPUT_NAME "molassembler"
  LIBRARY_OUTPUT_DIRECTORY "molassembler"
)
target_link_libraries(molassembler_python PRIVATE molassembler Eigen3::Eigen)
install(CODE "execute_process(WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${PYTHON_EXECUTABLE} -m pip install --prefix=${CMAKE_INSTALL_PREFIX} --upgrade .)")

# Python tests
if(SCINE_BUILD_TESTS)
  add_test(
    NAME MASMPythonBindings
    COMMAND ${PYTHON_EXECUTABLE} -B -m pytest ${CMAKE_CURRENT_SOURCE_DIR} --junitxml=${CMAKE_CURRENT_BINARY_DIR}/pytest_report.xml
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()

# Add the setup file in order to enable distribution and installation using pip
file(
  COPY ${CMAKE_CURRENT_SOURCE_DIR}/setup.py
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)

file(
  WRITE ${CMAKE_CURRENT_BINARY_DIR}/molassembler/__init__.py
  "from .molassembler import *\ndel molassembler\n"
)

file(
  COPY ${CMAKE_CURRENT_SOURCE_DIR}/sphinx
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)

configure_file(
  ${CMAKE_CURRENT_BINARY_DIR}/sphinx/Makefile.in
  ${CMAKE_CURRENT_BINARY_DIR}/sphinx/Makefile
  @ONLY
)
