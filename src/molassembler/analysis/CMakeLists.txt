# Early return in case analysis binaries are not to be compiled
if(NOT MOLASSEMBLER_ANALYSIS_BINARIES)
  return()
endif()

if(NOT TARGET Boost::regex
  OR NOT TARGET Boost::program_options
  OR NOT TARGET Boost::system
  OR NOT TARGET Boost::filesystem
)
  find_package(Boost REQUIRED
    COMPONENTS regex program_options system filesystem
  )
endif()

set(analysisFiles
  ${CMAKE_CURRENT_SOURCE_DIR}/BenchmarkEigenvalues.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/BenchmarkGraphAlgorithms.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/BenchmarkIsomorphisms.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/BenchmarkRefinement.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/BenchmarkRefinementFunctions.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DirectedConformerGeneration.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Gor1Explainer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Interpret.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/MakeIsomorphism.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/RankingExplainer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/TraceRefinement.cpp
)

foreach(filename ${analysisFiles})
  get_filename_component(targetName ${filename} NAME_WE)

  # add the executable and specify compilation and linking options
  add_executable(${targetName} ${filename})

  target_link_libraries(${targetName}
    PRIVATE
      Boost::boost
      Boost::program_options
      Boost::regex
      Boost::system
      Boost::filesystem
      molassembler
      chemical_symmetries
      gor1
  )

  target_compile_definitions(${targetName} PRIVATE BOOST_ALL_DYN_LINK)
  target_compile_options(${targetName} PRIVATE ${MOLASSEMBLER_CXX_FLAGS})
endforeach()
