# Early return in case analysis binaries are not to be compiled
if(NOT MOLASSEMBLER_ANALYSIS_BINARIES)
  return()
endif()

project(molassemblerAnalyses)

find_package(Boost REQUIRED COMPONENTS regex program_options system)

set(analysisFiles
  ${CMAKE_CURRENT_SOURCE_DIR}/BenchmarkEigenvalues.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/BenchmarkGraphAlgorithms.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Editor.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Gor1Explainer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/MakeIsomorphism.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/RankingExplainer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/RaytraceRefinement.cpp
)

foreach(filename ${analysisFiles})
  set(targetName "")
  get_filename_component(targetName ${filename} NAME_WE)
  set(targetName "analysis_${targetName}")

  # add the executable and specify compilation and linking options
  add_executable(${targetName} ${filename})

  target_link_libraries(${targetName}
    PRIVATE
      Boost::boost
      Boost::program_options
      Boost::regex
      Boost::system
      molassembler
      gor1
  )

  target_compile_options(${targetName} PRIVATE ${MOLASSEMBLER_CXX_FLAGS})
endforeach()
