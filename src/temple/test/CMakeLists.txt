cmake_minimum_required(VERSION 3.0.2)

project(templeTestProject)

set(Boost_USE_STATIC_LIBS OFF)
find_package(Boost REQUIRED COMPONENTS unit_test_framework)
find_package(Eigen3 3.1.1 REQUIRED)

file(GLOB testFiles "*.cpp")
foreach(filename ${testFiles})
  # get a target name from the filename
  set(targetName "")
  get_filename_component(targetName ${filename} NAME_WE)
  set(targetName "templeTest_${targetName}")

  # add the executable and specify compilation and linking options
  add_executable(${targetName} ${filename})

  target_compile_options(${targetName} PUBLIC -Wall -Wpedantic -Wextra)
  if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    target_compile_options(${targetName} PUBLIC -O3 -march=native)
  else()
    target_compile_options(${targetName} PUBLIC -g)
    # target_compile_options(
    #  ${targetName} PUBLIC -g -fsanitize=address -fsanitize=undefined
    #)
    # set_target_properties(
    #  ${targetName} PROPERTIES
    #    LINK_FLAGS "-fsanitize=address -fsanitize=undefined"
    #)
  endif()

  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    target_compile_options(
      ${targetName}
      PRIVATE -fconstexpr-steps=1999999999 -fconstexpr-backtrace-limit=0
    )
  endif()

  target_include_directories(${targetName} PRIVATE ${Boost_INCLUDE_DIRS})
  target_link_libraries(
    ${targetName}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    temple
  )

  # add target to tests
  add_test(NAME ${targetName} COMMAND ${targetName})

  # set folder property
  set_property(TARGET ${targetName} PROPERTY FOLDER "Tests")
endforeach()
