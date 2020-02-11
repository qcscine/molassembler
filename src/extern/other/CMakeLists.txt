install(
  FILES eigen.txt boost.txt
  DESTINATION share/doc/molassembler/licenses
)

if(SCINE_BUILD_PYTHON_BINDINGS)
  install(
    FILES pybind.txt
    DESTINATION share/doc/molassembler/licenses
  )
endif()
