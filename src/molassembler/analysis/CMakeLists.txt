project(molassemblerAnalyses)

find_package(Boost REQUIRED COMPONENTS program_options system)

# Analysis targets. Build only if Debug
if(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  file(GLOB analysisFiles "*.cpp")
  foreach(filename ${analysisFiles})
    set(targetName "")
    get_filename_component(targetName ${filename} NAME_WE)
    set(targetName "analysis_${targetName}")

    # add the executable and specify compilation and linking options
    add_executable(${targetName} ${filename})
    #set_target_properties(
    #  ${targetName} PROPERTIES
    #  COMPILE_FLAGS ${WARNING_FLAGS}
    #  LINK_FLAGS "-fsanitize=address -fsanitize=undefined"
    #)
    target_compile_options(
      ${targetName}
        PRIVATE -Wall -Wpedantic -Wextra -Og -g
    )

    target_include_directories(${targetName} PRIVATE ${Boost_INCLUDE_DIRS})
    target_link_libraries(
      ${targetName}
      ${Boost_PROGRAM_OPTIONS_LIBRARY}
      ${Boost_SYSTEM_LIBRARY}
      molassemblerStatic
    )
  endforeach()
else()
  # Build the benchmark anyway
  add_executable(analysis_BenchmarkGraphAlgorithms BenchmarkGraphAlgorithms.cpp)
  target_compile_options(
    analysis_BenchmarkGraphAlgorithms
      PUBLIC -O3 -g -march=native -Wall -Wpedantic -Wextra
  )
  target_link_libraries(
    analysis_BenchmarkGraphAlgorithms
      ${Boost_PROGRAM_OPTIONS_LIBRARY}
      ${Boost_SYSTEM_LIBRARY}
      molassemblerStatic
  )
endif()
