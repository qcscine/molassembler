project(molassemblerTests)

find_package(Boost REQUIRED COMPONENTS unit_test_framework regex)

# Testing targets
# get all test cpp files
file(GLOB testFiles "*.cpp")
foreach(filename ${testFiles})
  # get a target name from the filename
  set(baseName "")
  get_filename_component(baseName ${filename} NAME_WE)

  set(targetName "test_${baseName}")

  # add the executable and specify compilation and linking options
  add_executable(${targetName} ${filename})

  if("${baseName}" STREQUAL "CycleData")
    message(STATUS "Adding sanitizers to ${targetName}")
    target_compile_options(${targetName}
      PUBLIC -g -fsanitize=address -fsanitize=undefined
    )
    set_target_properties(${targetName} PROPERTIES
      LINK_FLAGS "-fsanitize=address -fsanitize=undefined"
    )
  endif()

  target_link_libraries(${targetName}
    PRIVATE
      Boost::boost
      Boost::regex
      Boost::unit_test_framework
      molassembler
  )

  # add target to tests
  add_test(NAME ${targetName} COMMAND ${targetName})
endforeach()

# Copy testing files
file(
  COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_files
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)
