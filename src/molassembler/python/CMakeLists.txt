include(ImportPybind11)
import_pybind11()
set(PYBIND11_PYTHON_VERSION ${PYTHONVERSION})

file(GLOB_RECURSE MOLASSEMBLER_PYTHON_CPPS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# Python module
pybind11_add_module(molassembler_python ${MOLASSEMBLER_PYTHON_CPPS})
set_target_properties(molassembler_python PROPERTIES
  SUFFIX ".so"
  LIBRARY_OUTPUT_NAME "scine_molassembler"
  LIBRARY_OUTPUT_DIRECTORY "scine_molassembler"
)
target_include_directories(molassembler_python SYSTEM PRIVATE ${PYBIND11_INCLUDE_DIR})
target_compile_options(molassembler_python PRIVATE ${MOLASSEMBLER_CXX_FLAGS})
target_link_libraries(molassembler_python PRIVATE ${MOLASSEMBLER_STATIC_TARGET} Eigen3::Eigen)

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
    WRITE ${CMAKE_CURRENT_BINARY_DIR}/scine_molassembler/__init__.py
    "from .scine_molassembler import *\n"
  )

  install(CODE
    "execute_process(
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMAND ${PYTHON_EXECUTABLE} -m pip install --prefix=${CMAKE_INSTALL_PREFIX} --upgrade .
      RESULT_VARIABLE retcode
    )
    if(NOT \$\{retcode\} EQUAL 0)
      message(FATAL_ERROR \"Non-zero exit code when installing Python module using pip\")
    endif()"
  )
else()
  message(WARNING "The Python bindings will not be installed since the python module pip is unavailable.")
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
      COMMENT "Generate python bindings documentation with sphinx"
    )

    install(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/
      DESTINATION share/doc/molassembler/python
      MESSAGE_NEVER
      OPTIONAL
    )

    find_python_module(doctest)
    if(PY_DOCTEST)
      find_python_module(scine_utilities)
      if(PY_SCINE_UTILITIES OR SCINE_UTILS_PYTHON_BINARY_PATH)
        add_test(
          NAME PythonBindingsDoctest
          COMMAND ${PYTHON_EXECUTABLE} -c "import sys, doctest, scine_molassembler; (f, _) = doctest.testmod(scine_molassembler.scine_molassembler); sys.exit(f > 0)"
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        )

        if(SCINE_UTILS_PYTHON_BINARY_PATH)
          set_tests_properties(PythonBindingsDoctest PROPERTIES
            ENVIRONMENT PYTHONPATH=${SCINE_UTILS_PYTHON_BINARY_PATH}:$ENV{PYTHONPATH}
          )
        endif()
      else()
        message(STATUS "Python documentation testing is disabled as the scine_utilities module cannot be found and is not being built in-tree")
      endif()
    else()
      message(STATUS "Python documentation testing is disabled as doctest cannot be found")
    endif()
  else()
    message(STATUS "The Python binding documentation will not be built as sphinx cannot be found")
  endif()
endif()
