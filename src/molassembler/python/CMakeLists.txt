include(ImportPybind11)
import_pybind11()
set(PYBIND11_PYTHON_VERSION ${PYTHONVERSION})

file(GLOB_RECURSE MOLASSEMBLER_PYTHON_CPPS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
file(GLOB_RECURSE MOLASSEMBLER_PYTHON_PYS ${CMAKE_CURRENT_SOURCE_DIR}/*.py)

# Python module
pybind11_add_module(molassembler_python ${MOLASSEMBLER_PYTHON_CPPS})
set_target_properties(molassembler_python PROPERTIES
  SUFFIX ".so"
  LIBRARY_OUTPUT_NAME "molassembler"
)
target_link_libraries(molassembler_python PRIVATE molassembler Eigen3::Eigen)

# Python tests
if(SCINE_BUILD_TESTS)
  add_test(
    NAME MASMPythonBindings
    COMMAND ${PYTHON_EXECUTABLE} -B -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/Tests --junitxml=${CMAKE_CURRENT_BINARY_DIR}/pytest_report.xml
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()
