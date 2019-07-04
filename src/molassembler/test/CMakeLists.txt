if(NOT TARGET Boost::unit_test_framework
    OR NOT TARGET Boost::regex
    OR NOT TARGET Boost::filesystem
)
  find_package(Boost REQUIRED COMPONENTS unit_test_framework regex filesystem)
endif()

set(MainTestFiles
  ${CMAKE_CURRENT_SOURCE_DIR}/AtomStereopermutator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/BondStereopermutator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Canonicalization.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CommonTrig.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Cycles.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DelibHelpers.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Descriptors.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Editing.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/InnerGraph.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Molecule.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/OrderDiscoveryHelper.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/RankingInformation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/RankingSpotCheck.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Serialization.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/VSEPR.cpp
)

set(DGTestFiles
  ${CMAKE_CURRENT_SOURCE_DIR}/DGDirectedConformerGenerator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DGExplicitGraph.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DGFixedPositions.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DGGraphShortestDistances.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DGImplicitGraph.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DGMetricMatrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DGRefinementProblem.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DGmain.cpp
)

# Testing targets
add_executable(MolassemblerTests ${MainTestFiles})
target_link_libraries(MolassemblerTests
  PRIVATE
    Boost::boost
    Boost::regex
    gor1
    molassembler
    stereopermutation
    Boost::unit_test_framework
    Boost::filesystem
)
target_compile_definitions(MolassemblerTests PRIVATE BOOST_ALL_DYN_LINK)
target_compile_options(MolassemblerTests PRIVATE ${MOLASSEMBLER_CXX_FLAGS})
add_test(
  NAME MASMMainTests
  COMMAND MolassemblerTests
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/../data/"
)

add_executable(DistanceGeometryTests ${DGTestFiles})
target_link_libraries(DistanceGeometryTests
  PRIVATE
    Boost::boost
    Boost::regex
    gor1
    molassembler
    stereopermutation
    Boost::unit_test_framework
    Boost::filesystem
)
target_compile_options(DistanceGeometryTests PRIVATE ${MOLASSEMBLER_CXX_FLAGS})
target_compile_definitions(DistanceGeometryTests PRIVATE BOOST_ALL_DYN_LINK)
add_test(
  NAME MASMDistanceGeometryTests
  COMMAND DistanceGeometryTests
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/../data/"
)
set_property(TARGET DistanceGeometryTests PROPERTY FOLDER "Tests")
