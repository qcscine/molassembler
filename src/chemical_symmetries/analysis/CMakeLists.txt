if(NOT MOLASSEMBLER_ANALYSIS_BINARIES)
  return()
endif()

project(ChemicalSymmetriesAnalysisBinaries)

find_package(Boost REQUIRED
  COMPONENTS unit_test_framework program_options system
)
find_package(Eigen3 3.1.1 REQUIRED)

set(analysisFiles
  ${CMAKE_CURRENT_SOURCE_DIR}/writeMappings.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/writeTransitionDotFiles.cpp
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
      Boost::system
      chemical_symmetries
  )

  target_compile_options(${targetName} PRIVATE ${MOLASSEMBLER_CXX_FLAGS})

  set_property(${targetName} writeMappings PROPERTY FOLDER "Analysis")
endforeach()
