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

  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    target_compile_options(
      ${targetName}
      PRIVATE -fconstexpr-steps=1999999999 -fconstexpr-backtrace-limit=0
    )
  endif()

  target_link_libraries(${targetName}
    PRIVATE
      Boost::boost
      Boost::unit_test_framework
      temple
  )

  # add target to tests
  add_test(NAME ${targetName} COMMAND ${targetName})

  # set folder property
  set_property(TARGET ${targetName} PROPERTY FOLDER "Tests")
endforeach()
