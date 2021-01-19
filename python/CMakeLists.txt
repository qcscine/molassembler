include(ImportPybind11)
import_pybind11()
set(PYBIND11_PYTHON_VERSION ${PYTHONVERSION})

file(GLOB_RECURSE MOLASSEMBLER_PYTHON_CPPS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# Python module
pybind11_add_module(scine_molassembler ${MOLASSEMBLER_PYTHON_CPPS})
set_target_properties(scine_molassembler PROPERTIES
  LIBRARY_OUTPUT_NAME "scine_molassembler"
  LIBRARY_OUTPUT_DIRECTORY "scine_molassembler"
)
target_include_directories(scine_molassembler SYSTEM PRIVATE ${pybind11_INCLUDE_DIR} ${PYBIND11_INCLUDE_DIR})
target_compile_options(scine_molassembler PRIVATE ${MOLASSEMBLER_CXX_FLAGS})
target_link_libraries(scine_molassembler PRIVATE ${MOLASSEMBLER_STATIC_TARGET} Eigen3::Eigen)

include(FindPythonModule)
find_python_module(pytest)
if(NOT PY_PYTEST)
  cmessage(WARNING "Python module will not be tested since the python module 'pytest' could not be found.")
endif()

# Python tests
if(SCINE_BUILD_TESTS AND PY_PYTEST)
  add_test(
    NAME MASMPythonBindings
    COMMAND ${PYTHON_EXECUTABLE} -B -m pytest ${CMAKE_CURRENT_SOURCE_DIR} --junitxml=${CMAKE_CURRENT_BINARY_DIR}/pytest_report.xml
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()

# Installation of the python module with pip
find_python_module(pip)
if(PY_PIP)
  file(
    COPY ${CMAKE_CURRENT_SOURCE_DIR}/pkginit.py
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scine_molassembler
  )
  file(
    RENAME ${CMAKE_CURRENT_BINARY_DIR}/scine_molassembler/pkginit.py
    ${CMAKE_CURRENT_BINARY_DIR}/scine_molassembler/__init__.py
  )

  # Copy utils library into the package if utils is shared. If utils
  # is static, it is subsumed and we don't need to copy it (and since we don't
  # use core, that doesn't get dragged along either). molassembler is always
  # linked statically into the python module, so no need to worry about it
  include(TargetLibName)
  target_lib_type(Scine::UtilsOS _utils_libtype)
  if("${_utils_libtype}" STREQUAL "SHARED_LIBRARY")
    # Overwrite the python module's RPATH to ensure it looks in its directory
    if(APPLE)
      set_target_properties(scine_molassembler PROPERTIES
        BUILD_WITH_INSTALL_RPATH ON
        INSTALL_RPATH "@loader_path;@loader_path/../lib"
      )
    elseif(UNIX)
      set_target_properties(scine_molassembler PROPERTIES
        BUILD_WITH_INSTALL_RPATH ON
        INSTALL_RPATH "\$ORIGIN;\$ORIGIN/../lib"
      )
    endif()

    # Copy Utils and Core (this works regardless of whether they are imported
    # or build in-tree)
    add_custom_command(TARGET scine_molassembler POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:Scine::UtilsOS> ${CMAKE_CURRENT_BINARY_DIR}/scine_molassembler
      COMMENT "Copying dependent shared libraries into python package directory"
    )
    target_lib_filename(Scine::UtilsOS _utils_name)
    set(molassembler_PY_DEPS ", \"${_utils_name}\"")
    unset(_utils_name)
    unset(_core_name)
  endif()
  unset(_utils_libtype)

  # Add setuptools file
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in
    ${CMAKE_CURRENT_BINARY_DIR}/setup.py
    @ONLY
  )
  unset(molassembler_PY_DEPS)

  install(CODE
    "execute_process(
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMAND ${PYTHON_EXECUTABLE} -m pip install --prefix=${CMAKE_INSTALL_PREFIX} --upgrade --no-deps .
      RESULT_VARIABLE retcode
    )
    if(NOT \$\{retcode\} EQUAL 0)
      message(FATAL_ERROR \"Non-zero exit code when installing Python module using pip\")
    endif()"
  )
else()
  cmessage(WARNING "The Python bindings will not be installed since the python module pip is unavailable.")
endif()

include(SphinxDocumentation)
scine_sphinx_documentation(
  TARGET scine_molassembler
  CONFIGURATION ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/conf.py.in
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/sphinx
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/doc
  LINK UtilsPythonDoc
  DOCTEST
  DOCTEST_REQUIRES scine_utilities
)
