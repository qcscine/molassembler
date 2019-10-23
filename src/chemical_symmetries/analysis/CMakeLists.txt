if(NOT MOLASSEMBLER_ANALYSIS_BINARIES)
  return()
endif()

if(NOT TARGET Boost::unit_test_framework OR NOT TARGET Boost::program_options OR NOT TARGET Boost::system)
  find_package(Boost REQUIRED
    COMPONENTS unit_test_framework program_options system
  )
endif()
find_package(Eigen3 3.1.1 REQUIRED)

if(MOLASSEMBLER_PARALLELIZE AND NOT TARGET OpenMP::OpenMP_CXX)
  find_package(OpenMP)
endif()

set(analysisFiles
  ${CMAKE_CURRENT_SOURCE_DIR}/thomson.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/csm.cpp
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
      $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
  )

  target_compile_options(${targetName} PRIVATE ${MOLASSEMBLER_CXX_FLAGS})
  target_compile_definitions(${targetName} PRIVATE BOOST_ALL_DYN_LINK)

  set_target_properties(${targetName} PROPERTIES FOLDER "Analysis")
endforeach()
