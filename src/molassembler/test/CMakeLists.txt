project(molassemblerTests)

set(Boost_USE_STATIC_LIBS OFF)
find_package(Boost REQUIRED COMPONENTS unit_test_framework regex)

# Testing targets
# get all test cpp files
file(GLOB testFiles "*.cpp")
foreach(filename ${testFiles})
  # get a target name from the filename
  set(targetName "")
  get_filename_component(targetName ${filename} NAME_WE)
  set(targetName "test_${targetName}")

  # add the executable and specify compilation and linking options
  add_executable(${targetName} ${filename})

  target_compile_options(${targetName} PRIVATE -Wall -Wpedantic -Wextra)
  if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    target_compile_options(${targetName} PRIVATE -O2 -march=native)
  else()
    target_compile_options(${targetName} PRIVATE -g)
    # target_compile_options(
    #  ${targetName} PUBLIC -g -fsanitize=address -fsanitize=undefined
    #)
    # set_target_properties(
    #  ${targetName} PROPERTIES
    #    LINK_FLAGS "-fsanitize=address -fsanitize=undefined"
    #)
  endif()

  target_link_libraries(${targetName}
    PRIVATE
      Boost::boost
      Boost::unit_test_framework
      Boost::regex
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
