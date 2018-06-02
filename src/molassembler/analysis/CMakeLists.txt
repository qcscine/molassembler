project(molassemblerAnalyses)

find_package(Boost REQUIRED COMPONENTS regex program_options system)

# Analysis targets. Build only if Debug
if(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  file(GLOB analysisFiles "*.cpp")
  foreach(filename ${analysisFiles})
    set(targetName "")
    get_filename_component(targetName ${filename} NAME_WE)
    set(targetName "analysis_${targetName}")

    # add the executable and specify compilation and linking options
    add_executable(${targetName} ${filename})

    #if("${targetName}" STREQUAL "analysis_RaytraceRefinement")
    #  message(STATUS "Adding sanitizers to ${targetName}")
    #  target_compile_options(${targetName}
    #    PUBLIC -g -fsanitize=address -fsanitize=undefined
    #  )
    #  set_target_properties(${targetName} PROPERTIES
    #    LINK_FLAGS "-fsanitize=address -fsanitize=undefined"
    #  )
    #endif()
    target_compile_options(${targetName} PRIVATE -Og -g)

    target_link_libraries(${targetName}
      PRIVATE
        Boost::boost
        Boost::program_options
        Boost::regex
        Boost::system
        molassembler
    )
  endforeach()
else()
  # Build the benchmark anyway
  add_executable(analysis_BenchmarkGraphAlgorithms BenchmarkGraphAlgorithms.cpp)
  target_compile_options(analysis_BenchmarkGraphAlgorithms
    PRIVATE -O3 -g -march=native
  )
  target_link_libraries(analysis_BenchmarkGraphAlgorithms
    PRIVATE Boost::program_options molassembler
  )
endif()
