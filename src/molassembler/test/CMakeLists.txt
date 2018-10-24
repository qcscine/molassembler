project(molassemblerTests)

if(NOT TARGET Boost::unit_test_framework
    OR NOT TARGET Boost::regex
    OR NOT TARGET Boost::filesystem
)
  find_package(Boost REQUIRED COMPONENTS unit_test_framework regex filesystem)
endif()

set(testFiles
  ${CMAKE_CURRENT_SOURCE_DIR}/BondStereopermutator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CommonTrig.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Cycles.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DGExplicitGraph.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DGGraphShortestDistances.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DGImplicitGraph.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DGMetricMatrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DGRefinementProblem.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DelibHelpers.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Descriptors.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Molecule.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/OrderDiscoveryHelper.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/PermutationState.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/RankingInformation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/RankingSpotCheck.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Serialization.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/StdlibTypeAlgorithms.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/VSEPR.cpp
)

# Testing targets
foreach(filename ${testFiles})
  # get a target name from the filename
  set(baseName "")
  get_filename_component(baseName ${filename} NAME_WE)

  set(targetName "test_${baseName}")

  # add the executable and specify compilation and linking options
  add_executable(${targetName} ${filename})

  target_link_libraries(${targetName}
    PRIVATE
      Boost::boost
      Boost::regex
      gor1
      molassembler
      stereopermutation
      Boost::unit_test_framework
      Boost::filesystem
  )

  target_compile_options(${targetName} PRIVATE ${MOLASSEMBLER_CXX_FLAGS})

  # add target to tests
  add_test(
    NAME ${targetName}
    COMMAND ${targetName}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/../data/"
  )
endforeach()
