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
target_include_directories(molassembler_python SYSTEM PRIVATE ${PYBIND11_INCLUDE_DIR})
target_compile_options(molassembler_python PRIVATE ${MOLASSEMBLER_CXX_FLAGS})
target_link_libraries(molassembler_python PRIVATE molassembler Eigen3::Eigen)

include(FindPythonModule)

find_python_module(pytest)
if(NOT PY_PYTEST)
  message(WARNING "Python module will not be tested since the python module pytest could not be found.")
endif()

# Python tests
if(SCINE_BUILD_TESTS AND PY_PYTEST)
  add_test(
    NAME MASMPythonBindings
    COMMAND ${PYTHON_EXECUTABLE} -B -m pytest ${CMAKE_CURRENT_SOURCE_DIR} --junitxml=${CMAKE_CURRENT_BINARY_DIR}/pytest_report.xml
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()

find_python_module(pip)
if(PY_PIP)
# Add the setup file in order to enable distribution and installation using pip
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in
    ${CMAKE_CURRENT_BINARY_DIR}/setup.py
    @ONLY
  )

  file(
    WRITE ${CMAKE_CURRENT_BINARY_DIR}/molassembler/__init__.py
    "from .molassembler import *\ndel molassembler\n"
  )

  install(CODE "execute_process(WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${PYTHON_EXECUTABLE} -m pip install --prefix=${CMAKE_INSTALL_PREFIX} --upgrade .)")
else()
  message(WARNING "Python will not be installed since the python module pip is unavailable.")
endif()

if(SCINE_BUILD_DOCS)
  find_python_module(sphinx)
  if(PY_SPHINX)
    file(
      COPY ${CMAKE_CURRENT_SOURCE_DIR}/sphinx
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
    )

    configure_file(
      ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/conf.py.in
      ${CMAKE_CURRENT_BINARY_DIR}/sphinx/conf.py
      @ONLY
    )

    file(GLOB_RECURSE MOLASSEMBLER_PYTHON_DOC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/*.rst)

    add_custom_target(PythonDocumentation
      ALL
      COMMAND ${PYTHON_EXECUTABLE} -m sphinx -b html -a -E "sphinx" "doc"
      DEPENDS molassembler_python
      SOURCES ${MOLASSEMBLER_PYTHON_DOC_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/conf.py
      COMMENT "Generate python documentation with sphinx"
    )

    install(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/html/
      DESTINATION share/molassembler/sphinx
      MESSAGE_NEVER
      OPTIONAL
    )
  endif()
endif()
